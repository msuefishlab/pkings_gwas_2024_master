#!/usr/bin/env python3
"""05_deserts.py
Call-desert mask (DESIGN.md section 3): per pair, a 100-kb bin is excluded when its callable SNPs that are
polymorphic in the pair, per callable site, fall below 10% of the genome-wide median for that pair.
Reads sfs/chunks/bins.*.tsv from `03_spectra.py --mode bins`; writes sfs/deserts.tsv (excluded bins) and
sfs/bin_density.tsv (all bins).
"""
import collections
import glob
import os
import re
import statistics
import sys

OUT = os.path.join(os.getcwd(), os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
SUFFIX = sys.argv[1] if len(sys.argv) > 1 else ""
pat = re.compile(r"bins\.\d+" + re.escape(SUFFIX) + r"\.tsv$")
v = collections.defaultdict(lambda: collections.Counter())
for path in [f for f in glob.glob(os.path.join(OUT, "sfs/chunks/bins.*.tsv")) if pat.search(f)]:
    with open(path) as fh:
        next(fh)
        for line in fh:
            pid, chrom, b, kind, n = line.rstrip("\n").split("\t")
            v[(pid, chrom, int(b))][kind] += int(n)

by_pair = collections.defaultdict(list)
for (pid, chrom, b), c in v.items():
    sites = c["callable_invariant"] + c["callable_snp_rows"]
    if sites:
        by_pair[pid].append((chrom, b, sites, c["polymorphic"], c["polymorphic"] / sites))
with open(os.path.join(OUT, f"sfs/deserts{SUFFIX}.tsv"), "w") as d, open(os.path.join(OUT, f"sfs/bin_density{SUFFIX}.tsv"), "w") as allf:
    d.write("pair_id\tchrom\tbin100k\tcallable_sites\tpolymorphic\tdensity\tmedian\n")
    allf.write("pair_id\tchrom\tbin100k\tcallable_sites\tpolymorphic\tdensity\n")
    for pid, rows in sorted(by_pair.items()):
        med = statistics.median(r[4] for r in rows)
        n_des = 0
        for chrom, b, sites, poly, dens in sorted(rows, key=lambda r: (int(r[0][3:]), r[1])):
            allf.write(f"{pid}\t{chrom}\t{b}\t{sites}\t{poly}\t{dens:.6f}\n")
            if dens < 0.10 * med:
                d.write(f"{pid}\t{chrom}\t{b}\t{sites}\t{poly}\t{dens:.6f}\t{med:.6f}\n")
                n_des += 1
        print(pid, "bins", len(rows), "median density", f"{med:.5f}", "deserts", n_des,
              "desert sites", sum(r[2] for r in rows if r[4] < 0.10 * med))
