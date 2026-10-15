#!/usr/bin/env python3
"""site_table_stream.py
Reduce one chunk of the all-sites VCF to (a) a table of SNP records and (b) compact summaries of
the invariant sites, which are far too many to write out one by one.

Reads, on stdin, `bcftools query` output with columns
    CHROM POS REF ALT  [AN_g AC_g for every group g in --groups]  [ExcHet_g for every --exchet group]
from a pass of `bcftools +fill-tags -- -S groups.tsv -t AN,AC,ExcHet`.

Writes
  stdout                           SNP rows: chrom pos ref alt aa dbin excl + the per-group columns
                                   (the caller bgzips it)
  <outdir>/invariant_hist.<id>.tsv per analysis pair: invariant-site counts by
                                   (100-kb bin, exon-distance bin, exclusion, outgroup callable,
                                   AN_bp, AN_tp), with AN capped at --an_cap. A projection to n <= cap
                                   haplotypes needs only "AN >= n", so the cap loses nothing.
  <outdir>/invariant_sums.<id>.tsv exact sums for the pixy check: per whole-population group
                                   sum AN(AN-1)/2 (pi comparisons) and per check pair sum AN_a*AN_b
                                   (dxy comparisons), the same sums over SNP records together with
                                   their differences (so pi and dxy close without re-reading the
                                   SNP table), and record counts
  <outdir>/index_snps.<id>.tsv     the SNP rows at the GWAS index SNPs, if any fall in this chunk

aa:   ancestral allele from the outgroup group OUT: "R" if OUT is fixed for REF, "A" if fixed for ALT,
      with AN_OUT >= --out_min_an; otherwise "." (unpolarized).
dbin: distance to the nearest exon: 0 inside, 1 <= 2 kb, 2 <= 5 kb, 3 <= 10 kb, 4 > 10 kb.
excl: 0 none, 1 GWAS extended peak +/-250 kb, 2 chr4 inversion-like polymorphism.
"""
import argparse
import collections
import sys

DBIN_EDGES = (2_000, 5_000, 10_000)


def load_bed(path, chrom, named=False):
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == chrom:
                out.append((int(f[1]), int(f[2]), f[3] if named and len(f) > 3 else ""))
    return sorted(out)


def to_int(x):
    return 0 if x in (".", "") else int(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--groups", required=True, help="comma-separated group order of the AN/AC columns")
    ap.add_argument("--exchet", required=True, help="comma-separated groups with an ExcHet column")
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--exclude", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--chunk_id", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--an_cap", type=int, default=60)
    ap.add_argument("--out_min_an", type=int, default=8)
    ap.add_argument("--index_snps", default="", help="peaks BED (col 6 = index SNP); their rows are also "
                    "written to <outdir>/index_snps.<id>.tsv")
    ap.add_argument("--bin", type=int, default=100_000,
                    help="bin size for invariant tallies (100 kb: fine enough for the call-desert mask, "
                         "and 1-Mb bootstrap blocks are unions of bins)")
    a = ap.parse_args()

    groups = a.groups.split(",")
    gi = {g: i for i, g in enumerate(groups)}
    exchet = a.exchet.split(",")
    out_i = gi["OUT"]

    pairs, checks = [], []
    with open(a.pairs) as fh:
        next(fh)
        for line in fh:
            pid, bp, tp, role = line.rstrip("\n").split("\t")
            (checks if role == "check" else pairs).append((pid, gi[bp], gi[tp]))
    whole = [g for g in groups if g.startswith("w")]

    index_pos = set()
    if a.index_snps:
        with open(a.index_snps) as fh:
            for line in fh:
                f = line.split()
                c = f[0] if f[0].startswith("chr") else "chr" + f[0]
                if c == a.chrom:
                    index_pos.add(int(f[5]))
    index_rows = []
    pi_snp_d = collections.Counter(); pi_snp_c = collections.Counter()
    dxy_snp_d = collections.Counter(); dxy_snp_c = collections.Counter()
    whole_i = [(g, gi[g]) for g in groups if g.startswith("w")]

    exons = load_bed(a.exons, a.chrom)
    excl = load_bed(a.exclude, a.chrom, named=True)
    ei = 0

    header = ["chrom", "pos", "ref", "alt", "aa", "dbin", "excl"]
    for g in groups:
        header += [f"AN_{g}", f"AC_{g}"]
    header += [f"ExcHet_{g}" for g in exchet]
    out = sys.stdout
    out.write("#" + "\t".join(header) + "\n")

    hist = collections.Counter()
    pi_comp = collections.Counter()
    dxy_comp = collections.Counter()
    cap = a.an_cap

    def expand(inv):
        """Fold one block's invariant-site patterns into the per-pair histograms and the sums."""
        for (block, dbin, xcode, rest), n in inv.items():
            cols = rest.split("\t")
            an = [to_int(cols[2 * i]) for i in range(len(groups))]
            out_ok = int(an[out_i] >= a.out_min_an)
            for g in whole:
                k = an[gi[g]]
                pi_comp[g] += n * k * (k - 1) // 2
            for pid, i, j in checks:
                dxy_comp[pid] += n * an[i] * an[j]
            for pid, i, j in pairs:
                hist[(pid, block, dbin, xcode, out_ok, min(an[i], cap), min(an[j], cap))] += n
        inv.clear()

    # invariant sites are tallied by their full AN pattern and expanded once per 100-kb bin,
    # which bounds memory to one bin's worth of distinct patterns
    inv = collections.Counter()          # (block, dbin, excl, rest-of-line) -> n sites
    cur_block = None
    n_snp = n_inv = n_skip = 0
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t", 4)
        pos = int(f[1])
        x = pos - 1                       # 0-based, to match the BED files
        # distance to the nearest exon (exons are merged and sorted; positions arrive in order)
        while ei < len(exons) and exons[ei][1] <= x:
            ei += 1
        if ei < len(exons) and exons[ei][0] <= x:
            dbin = 0
        else:
            d = min(x - exons[ei - 1][1] + 1 if ei > 0 else 10**9,
                    exons[ei][0] - x if ei < len(exons) else 10**9)
            dbin = 1 + sum(d > e for e in DBIN_EDGES)
        xcode = 0
        for s, e, name in excl:
            if s <= x < e:
                xcode = 2 if name == "chr4_inversion" else 1
                break
        ref, alt = f[2], f[3]
        block = x // a.bin
        if block != cur_block:
            expand(inv)
            cur_block = block
        if alt == ".":
            n_inv += 1
            inv[(block, dbin, xcode, f[4])] += 1
            continue
        if len(ref) != 1 or len(alt) != 1:
            n_skip += 1
            continue
        n_snp += 1
        cols = f[4].split("\t")
        an_o, ac_o = to_int(cols[2 * out_i]), to_int(cols[2 * out_i + 1])
        aa = "."
        if an_o >= a.out_min_an:
            aa = "R" if ac_o == 0 else ("A" if ac_o == an_o else ".")
        row = f"{f[0]}\t{pos}\t{ref}\t{alt}\t{aa}\t{dbin}\t{xcode}\t{f[4]}\n"
        out.write(row)
        if pos in index_pos:
            index_rows.append(row)
        # pi and dxy contributions of this SNP (pixy's definitions: differing pairs / pairs compared)
        for g, i in whole_i:
            an, ac = to_int(cols[2 * i]), to_int(cols[2 * i + 1])
            pi_snp_d[g] += ac * (an - ac)
            pi_snp_c[g] += an * (an - 1) // 2
        for pid, i, j in checks:
            ai, ci = to_int(cols[2 * i]), to_int(cols[2 * i + 1])
            aj, cj = to_int(cols[2 * j]), to_int(cols[2 * j + 1])
            dxy_snp_d[pid] += ci * (aj - cj) + cj * (ai - ci)
            dxy_snp_c[pid] += ai * aj

    expand(inv)

    with open(f"{a.outdir}/invariant_hist.{a.chunk_id}.tsv", "w") as fh:
        fh.write("pair_id\tchrom\tblock\tdbin\texcl\tout_ok\tan_bp\tan_tp\tn\n")
        for (pid, block, dbin, xcode, out_ok, i, j), n in sorted(hist.items()):
            fh.write(f"{pid}\t{a.chrom}\t{block}\t{dbin}\t{xcode}\t{out_ok}\t{i}\t{j}\t{n}\n")
    with open(f"{a.outdir}/invariant_sums.{a.chunk_id}.tsv", "w") as fh:
        fh.write("kind\tname\tchrom\tvalue\n")
        for g, v in sorted(pi_comp.items()):
            fh.write(f"pi_comp_invariant\t{g}\t{a.chrom}\t{v}\n")
        for pid, v in sorted(dxy_comp.items()):
            fh.write(f"dxy_comp_invariant\t{pid}\t{a.chrom}\t{v}\n")
        for g in sorted(pi_snp_c):
            fh.write(f"pi_diffs_snp\t{g}\t{a.chrom}\t{pi_snp_d[g]}\n")
            fh.write(f"pi_comp_snp\t{g}\t{a.chrom}\t{pi_snp_c[g]}\n")
        for pid in sorted(dxy_snp_c):
            fh.write(f"dxy_diffs_snp\t{pid}\t{a.chrom}\t{dxy_snp_d[pid]}\n")
            fh.write(f"dxy_comp_snp\t{pid}\t{a.chrom}\t{dxy_snp_c[pid]}\n")
        fh.write(f"n_records\tsnp\t{a.chrom}\t{n_snp}\n")
        fh.write(f"n_records\tinvariant\t{a.chrom}\t{n_inv}\n")
        fh.write(f"n_records\tskipped_non_snv\t{a.chrom}\t{n_skip}\n")

    if index_rows:
        with open(f"{a.outdir}/index_snps.{a.chunk_id}.tsv", "w") as fh:
            fh.write("#" + "\t".join(header) + "\n")
            fh.writelines(index_rows)


if __name__ == "__main__":
    main()
