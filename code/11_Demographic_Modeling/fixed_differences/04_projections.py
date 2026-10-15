#!/usr/bin/env python3
"""04_projections.py
Projection sizes per pair and side (DESIGN.md section 3): n = the 10th percentile of AN over masked,
polarized SNPs (so 90% of them have AN >= n), rounded down to an even number and capped at 60.
Reads sfs/chunks/an_hist.*.tsv from `03_spectra.py --mode an_hist`; writes sfs/projections.tsv.
Pass --cap 40 for the pre-specified run-time fallback.
"""
import argparse
import collections
import glob
import os
import re

OUT = os.path.join(os.getcwd(), os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
ap = argparse.ArgumentParser()
ap.add_argument("--cap", type=int, default=60)
ap.add_argument("--suffix", default="")
a = ap.parse_args()

h = collections.defaultdict(collections.Counter)
pat = re.compile(r"an_hist\.\d+" + re.escape(a.suffix) + r"\.tsv$")
files = sorted(f for f in glob.glob(os.path.join(OUT, "sfs/chunks/an_hist.*.tsv")) if pat.search(f))
for path in files:
    with open(path) as fh:
        next(fh)
        for line in fh:
            pid, side, an, n = line.rstrip("\n").split("\t")
            h[(pid, side)][int(an)] += int(n)


def p10(counter):
    total = sum(counter.values())
    acc = 0
    for an in sorted(counter):
        acc += counter[an]
        if acc > 0.10 * total:        # the smallest AN with more than 10% of SNPs at or below it
            return an, total
    return max(counter), total


rows = {}
for (pid, side), c in h.items():
    an, total = p10(c)
    n = min(a.cap, an - an % 2)
    kept = sum(v for k, v in c.items() if k >= n) / total
    rows.setdefault(pid, {})[side] = (n, an, kept, total)
with open(os.path.join(OUT, f"sfs/projections{a.suffix}.tsv"), "w") as o:
    o.write("pair_id\tn_bp\tn_tp\tp10_AN_bp\tp10_AN_tp\tkept_bp\tkept_tp\tn_snps\n")
    for pid in sorted(rows):
        b, t = rows[pid]["bp"], rows[pid]["tp"]
        o.write(f"{pid}\t{b[0]}\t{t[0]}\t{b[1]}\t{t[1]}\t{b[2]:.3f}\t{t[2]:.3f}\t{b[3]}\n")
        print(pid, "n_bp", b[0], "n_tp", t[0], f"kept {b[2]:.3f}/{t[2]:.3f}")
