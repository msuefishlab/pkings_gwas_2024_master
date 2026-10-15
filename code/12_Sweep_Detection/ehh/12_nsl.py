#!/usr/bin/env python3
"""
12_nsl.py
nSL (number of segregating sites by length) per lineage, one chromosome at a time.

nSL is a haplotype statistic like iHS, but it counts segregating sites rather than integrating
over physical distance, which makes it comparatively insensitive to variation in the
recombination rate. No recombination map exists for this assembly, so the haplotype statistics
cannot be calibrated against a map; nSL is included as the robustness check for that gap
(DESIGN.md section 4, secondary statistics).

Like iHS, nSL contrasts the derived and ancestral haplotypes at each site, so it has little
power where the derived allele is close to fixation. The applicability map in
evidence/applicability_map.tsv records where that applies.

Scores are standardized within derived-allele-count bins, then |nSL| is written: a sweep can
lengthen either the derived or the ancestral haplotype depending on polarization.

Usage (inside dadi_tools.sif):
  python3 12_nsl.py --vcf <per-chromosome polarized VCF> --chrom chr6 --out <out.tsv.gz>
"""

import argparse
import gzip
import sys

import numpy as np

try:
    import allel
except ImportError:
    sys.exit("scikit-allel is required (run this inside dadi_tools.sif)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="polarized, phased per-chromosome VCF")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out", required=True, help="output .tsv.gz (chrom, pos, value)")
    ap.add_argument("--min-maf", type=float, default=0.05)
    ap.add_argument("--bins", type=int, default=20,
                    help="derived-allele-count bins for standardization")
    args = ap.parse_args()

    print(f"Reading {args.vcf}")
    cb = allel.read_vcf(args.vcf, fields=["variants/POS", "variants/REF", "variants/ALT",
                                          "variants/AA", "calldata/GT"])
    if cb is None or "calldata/GT" not in cb:
        sys.exit(f"no genotypes read from {args.vcf}")

    pos = cb["variants/POS"]
    gt = allel.GenotypeArray(cb["calldata/GT"])
    h = gt.to_haplotypes()
    print(f"  {h.shape[0]} variants x {h.shape[1]} haplotypes")

    # Orient alleles so that 1 = derived, using INFO/AA where it is available. nSL is defined
    # on derived versus ancestral haplotypes, so an unpolarized site would contribute a score
    # with the wrong sign.
    aa = cb.get("variants/AA")
    ref = cb["variants/REF"]
    alt = cb["variants/ALT"][:, 0]
    if aa is not None:
        aa = np.array([str(a) if a is not None else "." for a in aa])
        flip = (aa == alt) & (aa != ".")
        keep_pol = (aa == ref) | (aa == alt)
        print(f"  ancestral allele known at {keep_pol.sum()} of {len(pos)} sites "
              f"({100.0 * keep_pol.sum() / max(len(pos), 1):.1f}%); {flip.sum()} flipped")
    else:
        flip = np.zeros(len(pos), dtype=bool)
        keep_pol = np.ones(len(pos), dtype=bool)
        print("  WARNING: no INFO/AA field; assuming REF is ancestral")

    hv = np.asarray(h)
    hv[flip] = 1 - hv[flip]

    ac = hv.sum(axis=1)
    n_hap = hv.shape[1]
    freq = ac / n_hap
    keep = keep_pol & (freq >= args.min_maf) & (freq <= 1 - args.min_maf)
    print(f"  {keep.sum()} sites after polarization and MAF >= {args.min_maf}")
    if keep.sum() < 100:
        sys.exit("too few sites to compute nSL")

    hk = allel.HaplotypeArray(hv[keep])
    pk = pos[keep]
    ack = ac[keep]

    print("  computing nSL")
    score = allel.nsl(hk, use_threads=False)

    # standardize within derived-allele-count bins, as for iHS
    std, _ = allel.standardize_by_allele_count(score, ack, n_bins=args.bins, diagnostics=False)
    val = np.abs(std)

    ok = np.isfinite(val)
    print(f"  {ok.sum()} finite standardized scores")
    with gzip.open(args.out, "wt") as out:
        out.write("chrom\tpos\tvalue\n")
        for p, v in zip(pk[ok], val[ok]):
            out.write(f"{args.chrom}\t{p}\t{v:.6f}\n")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
