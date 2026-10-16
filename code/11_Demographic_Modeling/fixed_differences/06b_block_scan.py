#!/usr/bin/env python3
"""06b_block_scan.py
Callable sites per 1-Mb block over the whole scan (every exon class and region, outside call deserts), per
pair and variant. E_hetNe weights each block by these lengths (DESIGN.md section 4: E_hetNe = sum_w L_w q(f_w),
with L_w the block's share of L_scan). 06_gather.py stores per-block lengths for the fitting mask only.

Reads sfs/chunks/spectra.<chunk><suffix>.pkl; writes sfs/<pair><suffix>.block_scan.tsv
(chrom, block, L_scan, L_scan_alt) and checks that the blocks sum to the pair's L_scan in sfs/<pair>.pkl.

Usage (inside dadi_tools.sif, from the repo root):
  python3 .../06b_block_scan.py [--suffix _x]
"""
import argparse
import collections
import glob
import os
import pickle
import re

OUT = os.path.join(os.getcwd(), os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
VARIANTS = (("", "spectra", "invariant"), ("_alt", "spectra_alt", "invariant_all"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    pat = re.compile(r"spectra\.\d+" + re.escape(a.suffix) + r"\.pkl$")
    files = sorted(f for f in glob.glob(os.path.join(OUT, "sfs/chunks/spectra.*.pkl")) if pat.search(f))
    L = collections.defaultdict(lambda: collections.Counter())        # (pid, var) -> {(chrom, blk): n}
    for path in files:
        with open(path, "rb") as fh:
            p = pickle.load(fh)
        chrom = p["chrom"]
        for var, spec_key, inv_key in VARIANTS:
            if spec_key not in p:            # runs made before the label-free spectra existed
                continue
            for pid, sp in p[spec_key].items():
                for (blk, c3, region), arr in sp.items():
                    L[(pid, var)][(chrom, blk)] += float(arr.sum())
            for (pid, blk, c3, reg), n in p[inv_key].items():
                L[(pid, var)][(chrom, blk)] += n
    for pid in sorted({k[0] for k in L}):
        with open(os.path.join(OUT, f"sfs/{pid}{a.suffix}.pkl"), "rb") as fh:
            pk = pickle.load(fh)
        for var, _, _ in VARIANTS:
            if f"L_scan{var}" not in pk:
                continue
            tot, want = sum(L[(pid, var)].values()), pk[f"L_scan{var}"]
            if abs(tot - want) > 1e-6 * want:
                raise SystemExit(f"{pid}{var}: blocks sum to {tot:.0f}, sfs pickle says L_scan {want:.0f}")
        keys = sorted(set(L[(pid, "")]) | set(L[(pid, "_alt")]), key=lambda k: (int(k[0][3:]), k[1]))
        with open(os.path.join(OUT, f"sfs/{pid}{a.suffix}.block_scan.tsv"), "w") as o:
            o.write("chrom\tblock\tL_scan\tL_scan_alt\n")
            for k in keys:
                o.write(f"{k[0]}\t{k[1]}\t{L[(pid, '')][k]:.0f}\t{L[(pid, '_alt')][k]:.0f}\n")
        print(pid, "blocks", len(keys), "L_scan", f"{sum(L[(pid, '')].values()):.0f}",
              "L_scan_alt", f"{sum(L[(pid, '_alt')].values()):.0f}")


if __name__ == "__main__":
    main()
