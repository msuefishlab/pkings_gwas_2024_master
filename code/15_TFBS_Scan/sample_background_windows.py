#!/usr/bin/env python3
"""
Sample random background genomic windows matched to GWAS peak loci.

For each peak locus, samples N non-overlapping windows of the same size from
the same chromosome, excluding the peak region itself and annotated CDS exons.
Outputs a manifest TSV used by run_background_tfbs.sh.

Usage:
  python sample_background_windows.py \\
    --fai input_data/00_Reference_Genome/<assembly>.fa.fai \\
    --gff input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.gff \\
    --n_windows 100 \\
    --seed 42 \\
    --output output_data/15_TFBS_Scan/background/background_windows_manifest.tsv

Peak loci are hardcoded based on the three analyzed GWAS peaks:
  sptbn4b : chr17:11642886-11742887  (100 kb)
  cntnap5 : chr16:10839999-11040000  (200 kb)
  itih3   : chr8:15736955-15836956   (100 kb)

Sample IDs are accepted via --samples (order matches locus order above), or
individually with --sample_sptbn4b / --sample_cntnap5 / --sample_itih3.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np


# ── Hardcoded peak definitions ─────────────────────────────────────────────────
PEAKS = {
    "sptbn4b": {"chr": "chr17", "start": 11642886, "end": 11742887},
    "cntnap5": {"chr": "chr16", "start": 10839999, "end": 11040000},
    "itih3":   {"chr": "chr8",  "start": 15736955, "end": 15836956},
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sample random background windows matched to GWAS peak loci."
    )
    p.add_argument("--fai", required=True, help="Genome FASTA index (.fai)")
    p.add_argument("--gff", required=True, help="Genome annotation GFF (for CDS exclusion)")
    p.add_argument(
        "--n_windows", type=int, default=100,
        help="Background windows to sample per locus (default: 100)",
    )
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p.add_argument(
        "--output", default="output_data/15_TFBS_Scan/background/background_windows_manifest.tsv",
        help="Output manifest TSV path",
    )
    p.add_argument(
        "--samples", nargs=3, metavar=("SPTBN4B_SAMPLE", "CNTNAP5_SAMPLE", "ITIH3_SAMPLE"),
        help="TP sample IDs for sptbn4b, cntnap5, itih3 (in that order)",
    )
    p.add_argument("--sample_sptbn4b", default=None, help="TP sample ID for sptbn4b locus")
    p.add_argument("--sample_cntnap5", default=None, help="TP sample ID for cntnap5 locus")
    p.add_argument("--sample_itih3",   default=None, help="TP sample ID for itih3 locus")
    p.add_argument(
        "--no_cds_exclusion", action="store_true",
        help="Skip CDS exclusion (useful if GFF is unavailable or slow to parse)",
    )
    return p.parse_args()


def resolve_sample_ids(args) -> dict:
    sample_map = {}
    loci = list(PEAKS.keys())
    if args.samples:
        for locus, sid in zip(loci, args.samples):
            sample_map[locus] = sid
    for locus in loci:
        attr_name = f"sample_{locus}"
        val = getattr(args, attr_name, None)
        if val:
            sample_map[locus] = val
    missing = [l for l in loci if l not in sample_map]
    if missing:
        print(
            f"[WARN] No sample ID specified for: {', '.join(missing)}. "
            "These loci will have sample_id=UNSET in the manifest. "
            "Edit the manifest before running run_background_tfbs.sh.",
            file=sys.stderr,
        )
        for locus in missing:
            sample_map[locus] = "UNSET"
    return sample_map


def read_fai(fai_path: str) -> dict:
    """Return {chrom: length} from a FASTA index file."""
    chrom_sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def read_cds_intervals(gff_path: str, chromosomes: set) -> dict:
    """
    Return {chrom: sorted list of (start, end)} for CDS features.
    Coordinates are 0-based half-open (converted from GFF 1-based inclusive).
    """
    cds_by_chrom: dict = {c: [] for c in chromosomes}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type not in ("CDS", "exon"):
                continue
            chrom = parts[0]
            if chrom not in cds_by_chrom:
                continue
            start = int(parts[3]) - 1  # GFF 1-based → 0-based
            end = int(parts[4])        # GFF inclusive → half-open
            cds_by_chrom[chrom].append((start, end))
    for chrom in cds_by_chrom:
        cds_by_chrom[chrom].sort()
    return cds_by_chrom


def overlaps_any(start: int, end: int, intervals: list) -> bool:
    """Check if [start, end) overlaps any interval in a sorted list."""
    for (ivs, ive) in intervals:
        if ivs >= end:
            break
        if ive > start:
            return True
    return False


def sample_windows_for_locus(
    locus: str,
    peak: dict,
    chrom_size: int,
    exclude_intervals: list,
    n_windows: int,
    rng: np.random.Generator,
    max_attempts: int = 100_000,
) -> list:
    """
    Sample n_windows non-overlapping windows of size (peak.end - peak.start)
    on peak.chr, excluding the peak itself and any interval in exclude_intervals.

    Returns list of (start, end) tuples (0-based half-open).
    """
    size = peak["end"] - peak["start"]
    chr_name = peak["chr"]

    # Build exclusion list: peak + CDS
    all_exclude = [(peak["start"], peak["end"])] + exclude_intervals
    all_exclude.sort()

    sampled = []
    attempts = 0
    while len(sampled) < n_windows and attempts < max_attempts:
        attempts += 1
        max_start = chrom_size - size
        if max_start <= 0:
            print(f"[WARN] Chromosome {chr_name} ({chrom_size} bp) is too small for a {size} bp window.", file=sys.stderr)
            break
        start = int(rng.integers(0, max_start))
        end = start + size

        if overlaps_any(start, end, all_exclude):
            continue

        # Also don't overlap previously sampled windows
        if overlaps_any(start, end, sampled):
            continue

        sampled.append((start, end))

    if len(sampled) < n_windows:
        print(
            f"[WARN] Only found {len(sampled)}/{n_windows} non-overlapping background windows "
            f"for {locus} on {chr_name} after {max_attempts} attempts.",
            file=sys.stderr,
        )
    return sampled


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    sample_map = resolve_sample_ids(args)

    # Read chromosome sizes
    chrom_sizes = read_fai(args.fai)
    print(f"[INFO] Read {len(chrom_sizes)} chromosomes from {args.fai}")

    # Identify chromosomes needed
    needed_chroms = {p["chr"] for p in PEAKS.values()}
    missing_chroms = needed_chroms - set(chrom_sizes)
    if missing_chroms:
        print(f"[ERROR] Chromosomes not found in FAI: {missing_chroms}", file=sys.stderr)
        sys.exit(1)

    # Read CDS intervals for relevant chromosomes
    if args.no_cds_exclusion:
        cds_by_chrom = {c: [] for c in needed_chroms}
        print("[INFO] Skipping CDS exclusion (--no_cds_exclusion)")
    else:
        print(f"[INFO] Reading CDS/exon features from {args.gff} ...")
        cds_by_chrom = read_cds_intervals(args.gff, needed_chroms)
        for c in needed_chroms:
            print(f"[INFO]   {c}: {len(cds_by_chrom[c])} CDS/exon intervals")

    # Prepare output
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for locus, peak in PEAKS.items():
        chrom = peak["chr"]
        chrom_size = chrom_sizes[chrom]
        sample_id = sample_map[locus]
        cds_intervals = cds_by_chrom.get(chrom, [])

        print(f"\n[INFO] Sampling {args.n_windows} background windows for {locus} ({chrom}, {peak['end']-peak['start']} bp) ...")
        windows = sample_windows_for_locus(
            locus=locus,
            peak=peak,
            chrom_size=chrom_size,
            exclude_intervals=cds_intervals,
            n_windows=args.n_windows,
            rng=rng,
        )

        for i, (start, end) in enumerate(windows):
            window_id = f"bg_{locus}_{i+1:04d}"
            rows.append({
                "window_id": window_id,
                "locus":     locus,
                "chr":       chrom,
                "start":     start,
                "end":       end,
                "size_bp":   end - start,
                "sample_id": sample_id,
                "region":    f"{chrom}:{start+1}-{end}",  # 1-based for bcftools/get_tfbs
            })

        print(f"[INFO]   Sampled {len(windows)} windows for {locus}")

    # Write manifest
    fieldnames = ["window_id", "locus", "chr", "start", "end", "size_bp", "sample_id", "region"]
    with open(out_path, "w", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n[INFO] Wrote {len(rows)} background windows to: {out_path}")
    print("[INFO] Edit sample_id column if needed before running run_background_tfbs.sh")


if __name__ == "__main__":
    main()
