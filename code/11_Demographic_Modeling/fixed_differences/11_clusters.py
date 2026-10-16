#!/usr/bin/env python3
"""11_clusters.py
Near-fixed differences, their 75-kb clusters, the Phase B gate and the parallel count (DESIGN.md sections 4,
8 and 10, Deviation 1), read straight from the site table one chunk at a time.

  --chunk N    reads site_table/snps.N.tsv.gz; writes parts/clusters.N.tsv (per pair, threshold and
               direction: the chunk's clusters, each with its SNP count and most differentiated SNP) and
               parts/tallies.N.tsv (delta counts and parallel-count tallies)
  --combine    merges the parts and writes, in FD_OUT,
                 fixed_diff_clusters.<pair>.bed      BP-derived, >= 0.95 in BP and <= 0.05 in TP
                 fixed_diff_clusters_0.9.<pair>.bed  BP-derived, >= 0.9 / <= 0.1 (Phase B fallback)
                 fixed_diff_cluster_counts.tsv       pair, threshold, direction: SNPs and clusters
                 phaseB_gate.tsv                     section 8 gate for pBAM/pTP1, pBP2/pTP2, pBIR/pTP2
                 parallel_count.tsv                  section 4 parallel count for Peaks 2, 3 and 6

SNPs, as 03_spectra.py includes them: biallelic SNVs of the site table; ExcHet >= 0.001 in every whole
population; AN >= the pair's projection size in both groups (sfs/projections<suffix>.tsv, the default run
first) and outside the pair's call deserts (sfs/deserts<suffix>.tsv); polarized by outgroup agreement, derived
count d = AC if aa is REF, AN - AC if ALT. A pair that no 03_spectra.py run projected falls back to
AN >= 80% of 2n in both groups (group_sizes.tsv) with no desert mask; --an_rule frac applies that AN rule to
every pair. Extended peaks +/-250 kb (region peak1..peak6) are left out; exons and the chr4 inversion stay
in (06_gather.py scope "nonpeak"). Pixy check pairs (role check) have no BP/TP sides and are skipped.

Near-fixed (exact integer tests, so a boundary value counts the same on both sides)
  BP_derived   d_bp/AN_bp >= t and d_tp/AN_tp <= 1 - t;  TP_derived the mirror image;  t = 0.95 or 0.9
  n_snps_delta the section 4 asymmetry count: delta >= t (BP_derived) or delta <= -t (TP_derived), with
               delta = f_bp - f_tp; this replaces the histogram count of fixed_diff_counts.tsv, which
               dropped delta = -0.95 exactly
Clusters: near-fixed SNPs no more than 75 kb from their neighbour join one cluster (single linkage). A chunk
writes its own clusters; --combine sorts them along each chromosome and joins neighbours closer than 75 kb,
which can only happen across a chunk boundary, so the result equals clustering each chromosome whole.
BED columns: chrom, start (first SNP, 0-based), end (last SNP), n_snps, centre_pos (1-based position of the
most differentiated SNP, largest delta, ties to the smallest position), max_delta (its delta).

Parallel count (section 4; Peaks 2 and 3 use BIK, BAVA and BIR, Peak 6 BIK and BAVA, all against TP2): SNPs
outside every extended peak where one allele is at least as frequent as at the index SNP in each BP population
and at most as frequent in TP2.
  label_free  every SNP, polarized or not; the allele is the one more frequent in BIK, at the index SNP too
              (ALT on a tie). The Peak 2 and Peak 6 index SNPs are unlabelled under Deviation 1 (Addendum 1b),
              so this is their only count.
  derived     polarized SNPs and the derived allele, where the index SNP is polarized (Peak 3)
Rows for whole populations (w*, headline) and the phenotype-pure lists (p*). AN and desert rules: each BP
population with TP2 as in its own pair (A2, A3, A4), so a counted SNP is callable in every pair involved.
Index-SNP counts come from site_table/index_snps.*.tsv at the positions in peaks.bed column 6.

Usage (inside dadi_tools.sif, from the repo root; 11_clusters.sb runs both steps):
  python3 .../11_clusters.py --chunk 17
  python3 .../11_clusters.py --combine
"""
import argparse
import collections
import glob
import gzip
import importlib.util
import math
import os
import re
from fractions import Fraction

ROOT = os.getcwd()
INP = os.path.join(ROOT, "input_data/11_Demographic_Modeling")
OUT = os.path.join(ROOT, os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
ST = os.path.join(OUT, "site_table")
PARTS = os.path.join(OUT, "parts")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("spectra", os.path.join(HERE, "03_spectra.py"))
spectra = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spectra)
to_int = spectra.to_int

MERGE_BP = 75_000
THRESHOLDS = (("0.95", 95), ("0.9", 90))       # (label, percent)
DIRECTIONS = ("BP_derived", "TP_derived")
GATE_MIN = 200
GATE_PAIRS = (("BP1", "pBAM", "pTP1"), ("BP2", "pBP2", "pTP2"), ("BP3", "pBIR", "pTP2"))
PARALLEL = {2: ("BIK", "BAVA", "BIR"), 3: ("BIK", "BAVA", "BIR"), 6: ("BIK", "BAVA")}
SAMPLES = (("whole", "w"), ("pure", "p"))


def ckey(c):
    m = re.match(r"chr(\d+)$", c)
    return (int(m.group(1)) if m else 10**6, c)


def load_pairs():
    pairs = []
    with open(os.path.join(INP, "pairs.tsv")) as fh:
        next(fh)
        for line in fh:
            pid, bp, tp, role = line.rstrip("\n").split("\t")
            if role != "check":
                pairs.append((pid, bp, tp))
    return pairs


def an_rules(pairs, mode, frac):
    """Per pair: (AN_min_bp, AN_min_tp, desert bins, description)."""
    sizes = {}
    with open(os.path.join(INP, "group_sizes.tsv")) as fh:
        next(fh)
        for line in fh:
            g, n = line.rstrip("\n").split("\t")[:2]
            sizes[g] = int(n)
    runs = []
    for path in glob.glob(os.path.join(OUT, "sfs/projections*.tsv")):
        sfx = os.path.basename(path)[len("projections"):-len(".tsv")]
        if os.path.exists(os.path.join(OUT, f"sfs/deserts{sfx}.tsv")):
            runs.append(sfx)
    proj = {}
    for sfx in sorted(runs, key=lambda s: (s != "", s)):        # default run first
        p, d = spectra.read_projections(sfx), spectra.read_deserts(sfx)
        for pid, n in p.items():
            proj.setdefault(pid, (n, d.get(pid, set()), sfx or "default"))
    rules = {}
    for pid, bp, tp in pairs:
        if pid in proj and mode == "projection":
            (n1, n2), des, run = proj[pid]
            rules[pid] = (n1, n2, des, f"projection:{run}")
        else:
            f = Fraction(str(frac))
            des, run = (proj[pid][1], proj[pid][2]) if pid in proj else (set(), None)
            rules[pid] = (math.ceil(f * 2 * sizes[bp]), math.ceil(f * 2 * sizes[tp]), des,
                          f"frac{frac}" + (f"+deserts:{run}" if run else ":no_deserts"))
    return rules, sizes


def group_rule(pairs, rules, sizes, bp, tp, frac):
    """AN minima and deserts for a BP group against tp: those of the pair (bp, tp), else the AN fraction."""
    for pid, b, t in pairs:
        if (b, t) == (bp, tp):
            return rules[pid][:3]
    f = Fraction(str(frac))
    return math.ceil(f * 2 * sizes[bp]), math.ceil(f * 2 * sizes[tp]), set()


def load_index():
    """{peak: site-table row at its index SNP} for the parallel-count peaks."""
    rows = {}
    for path in glob.glob(os.path.join(ST, "index_snps.*.tsv")):
        if path.endswith(".smoke.tsv"):
            continue
        with open(path) as fh:
            header = fh.readline().lstrip("#").rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(header, line.rstrip("\n").split("\t")))
                rows[(f["chrom"], int(f["pos"]))] = f
    out = {}
    with open(os.path.join(ROOT, "input_data/10_Introgression/peaks.bed")) as fh:
        for line in fh:
            f = line.split()
            c = f[0] if f[0].startswith("chr") else "chr" + f[0]
            pk = int(f[3])
            if pk in PARALLEL:
                out[pk] = (c, int(f[5]), rows.get((c, int(f[5]))))
    return out


def parallel_setup(pairs, rules, sizes, frac, ix=None):
    """One test per (peak, sample): groups, AN minima, deserts and the index SNP's oriented counts."""
    tests = []
    for pk, (chrom, pos, row) in sorted(load_index().items()):
        for sample, s in SAMPLES:
            bps, tp = [s + p for p in PARALLEL[pk]], s + "TP2"
            gr = [group_rule(pairs, rules, sizes, g, tp, frac) for g in bps]
            an_min = [r[0] for r in gr] + [max(r[1] for r in gr)]
            des = set().union(*(r[2] for r in gr))
            t = {"peak": pk, "sample": sample, "chrom": chrom, "pos": pos, "groups": bps + [tp],
                 "an_min": an_min, "des": des, "lf": None, "der": None, "aa": "missing"}
            if row is not None:
                cnt = [(to_int(row["AC_" + g]), to_int(row["AN_" + g])) for g in bps + [tp]]
                t["aa"] = row["aa"]
                if all(n for _, n in cnt):
                    alt = 2 * cnt[0][0] >= cnt[0][1]            # allele more frequent in BIK
                    t["lf"] = [(a if alt else n - a, n) for a, n in cnt]
                    if row["aa"] != ".":
                        t["der"] = [(a if row["aa"] == "R" else n - a, n) for a, n in cnt]
            if ix is not None:
                t["cols"] = [(ix["AN_" + g], ix["AC_" + g]) for g in bps + [tp]]
            tests.append(t)
    return tests


def rect(cnt, ref):
    """At least as frequent as ref in every BP group, at most as frequent in the TP group (last)."""
    for (a, n), (a0, n0) in zip(cnt[:-1], ref[:-1]):
        if a * n0 < a0 * n:
            return False
    (a, n), (a0, n0) = cnt[-1], ref[-1]
    return a * n0 <= a0 * n


def delta_of(db, nb, dt, nt, direction):
    d = Fraction(db, nb) - Fraction(dt, nt)
    return d if direction == "BP_derived" else -d


def run_chunk(a):
    pairs = load_pairs()
    rules, sizes = an_rules(pairs, a.an_rule, a.an_frac)
    open_cl, done = {}, []          # (pid, thr, direction) -> [first, last, n, centre, db, nb, dt, nt, delta]
    n_delta = collections.Counter()
    par = collections.defaultdict(lambda: [0, 0])     # (peak, sample, orientation) -> [hits, tested]
    regs, chrom, last_pos = None, None, -1

    def add(key, pos, db, nb, dt, nt):
        c = open_cl.get(key)
        d = delta_of(db, nb, dt, nt, key[2])
        if c is not None and pos - c[1] <= MERGE_BP:
            c[1] = pos
            c[2] += 1
            if d > c[8]:                   # ties keep the earlier (smaller) position
                c[3:9] = [pos, db, nb, dt, nt, d]
        else:
            if c is not None:
                done.append((key, c))
            open_cl[key] = [pos, pos, 1, pos, db, nb, dt, nt, d]

    with gzip.open(os.path.join(ST, f"snps.{a.chunk}.tsv.gz"), "rt") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        exchet_cols = [ix[c] for c in header if c.startswith("ExcHet_")]
        cols = [(pid, ix["AN_" + bp], ix["AC_" + bp], ix["AN_" + tp], ix["AC_" + tp], rules[pid])
                for pid, bp, tp in pairs]
        tests = parallel_setup(pairs, rules, sizes, a.an_frac, ix)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if regs is None:
                chrom = f[0]
                regs = spectra.load_regions(chrom)
            pos = int(f[1]); x = pos - 1
            if pos < last_pos:
                raise SystemExit(f"chunk {a.chunk}: positions not sorted at {chrom}:{pos}")
            last_pos = pos
            # 03_spectra.py: paralog filter, then region from the exclusion code
            if any(f[i] not in (".", "") and float(f[i]) < spectra.EXCHET_MIN for i in exchet_cols):
                continue
            xcode = int(f[6])
            region = "none" if xcode == 0 else spectra.region_of(regs, x)
            if region.startswith("peak"):
                continue
            aa = f[4]
            bkey = (chrom, x // 100_000)
            if aa != ".":
                for pid, i_an, i_ac, j_an, j_ac, (n1, n2, des, _) in cols:
                    nb, nt = to_int(f[i_an]), to_int(f[j_an])
                    if nb < n1 or nt < n2 or not nb or not nt or bkey in des:
                        continue
                    ab, at = to_int(f[i_ac]), to_int(f[j_ac])
                    db, dt = (ab, at) if aa == "R" else (nb - ab, nt - at)
                    if 10 * db < 9 * nb and 10 * dt < 9 * nt:
                        continue                  # not near-fixed at either threshold
                    for thr, t in THRESHOLDS:
                        if 100 * db >= t * nb and 100 * dt <= (100 - t) * nt:
                            direction, num = "BP_derived", db * nt - dt * nb
                        elif 100 * dt >= t * nt and 100 * db <= (100 - t) * nb:
                            direction, num = "TP_derived", dt * nb - db * nt
                        else:
                            continue
                        add((pid, thr, direction), pos, db, nb, dt, nt)
                        if 100 * num >= t * nb * nt:          # |delta| >= t on the same side
                            n_delta[(pid, thr, direction)] += 1
            for t in tests:
                if t["lf"] is None or bkey in t["des"]:
                    continue
                cnt = []
                for (ci, cj), m in zip(t["cols"], t["an_min"]):
                    n = to_int(f[ci])
                    if n < m or not n:
                        break
                    cnt.append((to_int(f[cj]), n))
                else:
                    k = (t["peak"], t["sample"])
                    par[k + ("label_free",)][1] += 1
                    if rect(cnt, t["lf"]) or rect([(n - c, n) for c, n in cnt], t["lf"]):
                        par[k + ("label_free",)][0] += 1
                    if t["der"] is not None and aa != ".":
                        par[k + ("derived",)][1] += 1
                        d = cnt if aa == "R" else [(n - c, n) for c, n in cnt]
                        if rect(d, t["der"]):
                            par[k + ("derived",)][0] += 1
    done.extend(open_cl.items())

    os.makedirs(PARTS, exist_ok=True)
    tmp = os.path.join(PARTS, f"clusters.{a.chunk}.tsv.tmp")
    with open(tmp, "w") as o:
        for pid, _, _ in pairs:
            n1, n2, _, desc = rules[pid]
            o.write(f"#rule\t{pid}\t{n1}\t{n2}\t{desc}\n")
        o.write("pair_id\tthreshold\tdirection\tchrom\tfirst_pos\tlast_pos\tn_snps\tcentre_pos\td_bp\tAN_bp\td_tp\tAN_tp\n")
        for (pid, thr, direction), c in sorted(done, key=lambda r: (r[0], r[1][0])):
            o.write(f"{pid}\t{thr}\t{direction}\t{chrom}\t" + "\t".join(str(v) for v in c[:8]) + "\n")
    os.replace(tmp, tmp[:-4])
    tmp = os.path.join(PARTS, f"tallies.{a.chunk}.tsv.tmp")
    with open(tmp, "w") as o:
        o.write("what\tkey1\tkey2\tkey3\tn\tn_tested\n")
        for (pid, thr, direction), n in sorted(n_delta.items()):
            o.write(f"delta\t{pid}\t{thr}\t{direction}\t{n}\t.\n")
        for (pk, sample, orient), (h, n) in sorted(par.items()):
            o.write(f"parallel\t{pk}\t{sample}\t{orient}\t{h}\t{n}\n")
    os.replace(tmp, tmp[:-4])
    print(f"chunk {a.chunk} ({chrom}): {len(done)} clusters, {sum(n_delta.values())} delta rows, "
          f"{sum(v[1] for k, v in par.items() if k[2] == 'label_free')} parallel tests")


def run_combine(a):
    pairs = load_pairs()
    rules, sizes = an_rules(pairs, a.an_rule, a.an_frac)
    n_chunks = sum(1 for _ in open(os.path.join(INP, "chunks.tsv"))) - 1
    files = {}
    for kind in ("clusters", "tallies"):
        pat = re.compile(kind + r"\.\d+\.tsv$")
        files[kind] = sorted(p for p in glob.glob(os.path.join(PARTS, f"{kind}.*.tsv")) if pat.search(p))
        if len(files[kind]) != n_chunks:
            raise SystemExit(f"{len(files[kind])} {kind} parts for {n_chunks} chunks; finish the --chunk tasks first")

    # every chunk must have used the same AN and desert rule, and the rule must still hold
    want = {pid: f"{r[0]}\t{r[1]}\t{r[3]}" for pid, r in rules.items()}
    rows = collections.defaultdict(list)
    for path in files["clusters"]:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == "#rule":
                    if want.get(f[1]) != "\t".join(f[2:5]):
                        raise SystemExit(f"{path}: rule for {f[1]} is {f[2:5]}, now {want.get(f[1])}; rerun the chunks")
                    continue
                if f[0] == "pair_id":
                    continue
                first, last, n, centre, db, nb, dt, nt = (int(v) for v in f[4:12])
                rows[(f[0], f[1], f[2])].append((f[3], first, last, n, centre, db, nb, dt, nt))
    n_delta = collections.Counter()
    par = collections.defaultdict(lambda: [0, 0])
    for path in files["tallies"]:
        with open(path) as fh:
            next(fh)
            for line in fh:
                what, k1, k2, k3, n, nt = line.rstrip("\n").split("\t")
                if what == "delta":
                    n_delta[(k1, k2, k3)] += int(n)
                else:
                    par[(int(k1), k2, k3)][0] += int(n)
                    par[(int(k1), k2, k3)][1] += int(nt)

    # join clusters across chunk boundaries
    merged = {}
    for key, rs in rows.items():
        out = []
        for chrom, first, last, n, centre, db, nb, dt, nt in sorted(rs, key=lambda r: (ckey(r[0]), r[1])):
            d = delta_of(db, nb, dt, nt, key[2])
            if out and out[-1][0] == chrom and first - out[-1][2] <= MERGE_BP:
                o = out[-1]
                o[2] = max(o[2], last)
                o[3] += n
                if (d, -centre) > (o[5], -o[4]):
                    o[4], o[5] = centre, d
            else:
                out.append([chrom, first, last, n, centre, d])
        merged[key] = out

    for pid, _, _ in pairs:
        for thr, name in (("0.95", f"fixed_diff_clusters.{pid}.bed"), ("0.9", f"fixed_diff_clusters_0.9.{pid}.bed")):
            with open(os.path.join(OUT, name), "w") as o:
                for chrom, first, last, n, centre, d in merged.get((pid, thr, "BP_derived"), []):
                    o.write(f"{chrom}\t{first - 1}\t{last}\t{n}\t{centre}\t{float(d):.4f}\n")

    with open(os.path.join(OUT, "fixed_diff_cluster_counts.tsv"), "w") as o:
        o.write("pair_id\tthreshold\tdirection\tn_snps\tn_clusters\tn_snps_delta\tan_rule\n")
        for pid, _, _ in pairs:
            n1, n2, _, desc = rules[pid]
            for thr, _ in THRESHOLDS:
                for direction in DIRECTIONS:
                    cl = merged.get((pid, thr, direction), [])
                    o.write(f"{pid}\t{thr}\t{direction}\t{sum(c[3] for c in cl)}\t{len(cl)}\t"
                            f"{n_delta[(pid, thr, direction)]}\t{desc} AN>={n1},{n2}\n")

    # section 8 gate: >= 200 clusters at 0.95/0.05, else at 0.9/0.1, else Phase B is not run for the lineage
    with open(os.path.join(OUT, "phaseB_gate.tsv"), "w") as o:
        o.write("lineage\tpair_id\tbp\ttp\tn_clusters_0.95\tn_clusters_0.9\tdecision\tan_rule\n")
        for lin, bp, tp in GATE_PAIRS:
            pid = next((p for p, b, t in pairs if (b, t) == (bp, tp)), None)
            if pid is None:
                raise SystemExit(f"no pair {bp}/{tp} in pairs.tsv")
            c95 = len(merged.get((pid, "0.95", "BP_derived"), []))
            c90 = len(merged.get((pid, "0.9", "BP_derived"), []))
            dec = "run_0.95" if c95 >= GATE_MIN else ("run_relaxed_0.9" if c90 >= GATE_MIN else "not_run")
            o.write(f"{lin}\t{pid}\t{bp}\t{tp}\t{c95}\t{c90}\t{dec}\t{rules[pid][3]}\n")
            print(f"Phase B {lin} {pid}: {c95} clusters at 0.95, {c90} at 0.9 -> {dec}")

    with open(os.path.join(OUT, "parallel_count.tsv"), "w") as o:
        o.write("peak\tsample\torientation\tchrom\tpos\tgroups\tidx_f_BIK\tidx_f_BAVA\tidx_f_BIR\tidx_f_TP2\t"
                "n_snps\tn_tested\tnote\n")
        for t in parallel_setup(pairs, rules, sizes, a.an_frac):
            for orient in ("label_free", "derived"):
                ref = t["lf"] if orient == "label_free" else t["der"]
                if ref is None:
                    fr = ["NA"] * 4
                    hits = tested = "NA"
                    note = {"missing": "index SNP not in the site table", ".": "index SNP unpolarized (Addendum 1b)"}.get(
                        t["aa"], "index SNP not called in every group")
                else:
                    fr = [f"{a_ / n:.4f}" for a_, n in ref]
                    if len(fr) == 3:
                        fr = fr[:2] + ["NA"] + fr[2:]
                    hits, tested = par[(t["peak"], t["sample"], orient)]
                    note = "allele more frequent in BIK" if orient == "label_free" else "derived allele"
                o.write(f"{t['peak']}\t{t['sample']}\t{orient}\t{t['chrom']}\t{t['pos']}\t{','.join(t['groups'])}\t"
                        + "\t".join(fr) + f"\t{hits}\t{tested}\t{note}\n")
    print(f"{len(pairs)} pairs, {n_chunks} chunks -> {OUT}")


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--chunk", help="chunk_id of site_table/snps.<chunk>.tsv.gz")
    g.add_argument("--combine", action="store_true")
    ap.add_argument("--an_rule", default="projection", choices=["projection", "frac"],
                    help="projection: 03_spectra.py's AN rule where a run projected the pair (else frac); "
                         "frac: AN >= --an_frac x 2n for every pair")
    ap.add_argument("--an_frac", type=float, default=0.8)
    a = ap.parse_args()
    if a.combine:
        run_combine(a)
    else:
        run_chunk(a)


if __name__ == "__main__":
    main()
