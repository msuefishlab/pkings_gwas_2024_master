#!/usr/bin/env python3
"""02_check_site_table.py
Gate before any statistic (DESIGN.md section 9.1). The site table must reproduce what pixy measured
on the same VCF, and the index SNPs must match the sweep analysis, before anything is built on it.

  1. pi per whole population, per chromosome and genome-wide, within 0.5% of pixy
     (input_data/11_Demographic_Modeling/pixy_reference_totals.tsv: pixy count_diffs and
     count_comparisons summed over windows)
  2. dxy for the six check pairs, same tolerance
  3. SNP record count, reported against the 51,073,851-site phased SNP set (not gated: the
     all-sites VCF also carries SNPs that the phasing inputs dropped)
  4. index-SNP ALT frequencies in the phenotype-pure lineage groups within 0.05 of
     output_data/12_Sweep_Detection/evidence/index_snp_facts.tsv, and the BP allele derived
     under the new outgroup calls

Writes output_data/11_Demographic_Modeling/fixed_differences/checks/site_table_checks.tsv and exits
non-zero if any gated check fails.

Usage (from the repo root, inside dadi_tools.sif or gwas_tools.sif):
  python3 code/11_Demographic_Modeling/fixed_differences/02_check_site_table.py
"""
import collections
import csv
import glob
import os
import sys

ROOT = os.getcwd()
ST = os.path.join(ROOT, os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"), "site_table")
OUT = os.path.join(ROOT, os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"), "checks")
TOL_REL = 0.005
TOL_FREQ = 0.05
PHASED_SNPS = 51_073_851
def to_int(x):
    return 0 if x in (".", "") else int(x)


LINEAGE_GROUP = {"BP1": "pBAM", "TP1": "pTP1", "BP2": "pBP2", "BP3": "pBIR", "TP2": "pTP2"}


def main():
    os.makedirs(OUT, exist_ok=True)
    chunks = sorted(glob.glob(os.path.join(ST, "invariant_sums.*.tsv")))
    chunks = [c for c in chunks if not c.endswith(".smoke.tsv")]
    n_expected = sum(1 for _ in open(os.path.join(ROOT, "input_data/11_Demographic_Modeling/chunks.tsv"))) - 1
    rows, failed = [], False

    def record(check, item, chrom, ours, ref, ok, note=""):
        nonlocal failed
        rows.append((check, item, chrom, ours, ref, "PASS" if ok else "FAIL", note))
        failed = failed or not ok

    record("chunks_present", "invariant_sums", "genome", len(chunks), n_expected, len(chunks) == n_expected)

    v = collections.defaultdict(int)          # (kind, name, chrom) -> value
    for path in chunks:
        with open(path) as fh:
            next(fh)
            for line in fh:
                kind, name, chrom, val = line.rstrip("\n").split("\t")
                v[(kind, name, chrom)] += int(val)
                v[(kind, name, "genome")] += int(val)

    ref = {}
    with open(os.path.join(ROOT, "input_data/11_Demographic_Modeling/pixy_reference_totals.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ref[(r["stat"], r["pop1"], r["pop2"], r["chrom"])] = (float(r["count_diffs"]), float(r["count_comparisons"]))
    chroms = sorted({k[2] for k in v if k[2] != "genome"}, key=lambda c: int(c[3:])) + ["genome"]

    # 1. pi
    pops = sorted({k[1][1:] for k in v if k[0] == "pi_comp_snp"})
    for p in pops:
        g = "w" + p
        for c in chroms:
            if ("pi", p, "", c) not in ref:
                continue
            d = v[("pi_diffs_snp", g, c)]
            n = v[("pi_comp_snp", g, c)] + v[("pi_comp_invariant", g, c)]
            rd, rn = ref[("pi", p, "", c)]
            ours, theirs = d / n, rd / rn
            rel = abs(ours / theirs - 1)
            record("pi", p, c, f"{ours:.7f}", f"{theirs:.7f}", rel <= TOL_REL,
                   f"rel {rel:.4f}; diffs ratio {d / rd:.4f}; comparisons ratio {n / rn:.4f}")

    # 2. dxy
    for k in sorted({k[1] for k in v if k[0] == "dxy_comp_snp"}):
        a, b = k.split("_")[1:3]
        p1, p2 = sorted((a, b))
        for c in chroms:
            if ("dxy", p1, p2, c) not in ref:
                continue
            d = v[("dxy_diffs_snp", k, c)]
            n = v[("dxy_comp_snp", k, c)] + v[("dxy_comp_invariant", k, c)]
            rd, rn = ref[("dxy", p1, p2, c)]
            ours, theirs = d / n, rd / rn
            rel = abs(ours / theirs - 1)
            record("dxy", f"{p1}-{p2}", c, f"{ours:.7f}", f"{theirs:.7f}", rel <= TOL_REL,
                   f"rel {rel:.4f}; diffs ratio {d / rd:.4f}; comparisons ratio {n / rn:.4f}")

    # 3. record counts (reported, not gated)
    n_snp = v[("n_records", "snp", "genome")]
    rows.append(("snp_records", "all-sites VCF", "genome", n_snp, PHASED_SNPS, "INFO",
                 f"excess over phased SNP set {n_snp - PHASED_SNPS}; invariant "
                 f"{v[('n_records', 'invariant', 'genome')]}; skipped non-SNV {v[('n_records', 'skipped_non_snv', 'genome')]}"))

    # 4. index SNPs
    idx = {}
    for path in glob.glob(os.path.join(ST, "index_snps.*.tsv")):
        if path.endswith(".smoke.tsv"):
            continue
        with open(path) as fh:
            header = fh.readline().lstrip("#").rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(header, line.rstrip("\n").split("\t")))
                idx[(f["chrom"], int(f["pos"]))] = f
    facts = os.path.join(ROOT, "output_data/12_Sweep_Detection/evidence/index_snp_facts.tsv")
    with open(facts) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["chrom"], int(r["pos"]))
            name = f"P{r['peak']} {r['lineage']}"
            if key not in idx:
                record("index_snp", name, r["chrom"], "missing", r["alt_freq"], False)
                continue
            f = idx[key]
            g = LINEAGE_GROUP[r["lineage"]]
            an, ac = to_int(f[f"AN_{g}"]), to_int(f[f"AC_{g}"])
            ours = ac / an if an else float("nan")
            ok = an > 0 and abs(ours - float(r["alt_freq"])) <= TOL_FREQ
            record("index_snp_freq", name, r["chrom"], f"{ours:.3f} (AN {an})", r["alt_freq"], ok)
            if r["role"] == "own":
                anc = {"R": f["ref"], "A": f["alt"]}.get(f["aa"], ".")
                ok = anc == r["ancestral_allele"] and anc != r["bp_allele"]
                record("index_snp_bp_derived", f"P{r['peak']}", r["chrom"], f"ancestral {anc}",
                       f"ancestral {r['ancestral_allele']}, BP {r['bp_allele']}", ok)

    out = os.path.join(OUT, "site_table_checks.tsv")
    with open(out, "w") as fh:
        fh.write("check\titem\tchrom\tours\treference\tresult\tnote\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    n_fail = sum(r[5] == "FAIL" for r in rows)
    for r in rows:
        if r[2] == "genome" or r[5] != "PASS":
            print("\t".join(str(x) for x in r))
    print(f"\n{len(rows)} checks, {n_fail} failed -> {out}")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
