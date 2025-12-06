#!/usr/bin/env python3
"""
Analyze BP vs TP TF motif differences at SPTBN4 (or any locus).

Input:
    A BED-like file produced by bedtools intersect -wa -wb, with columns:
        0: motif_chr
        1: motif_start (0-based)
        2: motif_end   (0-based, half-open)
        3: motif_id
        4: fimo_score
        5: haplotype    ("BP" or "TP")
        6: fimo_pval
        7: snp_chr
        8: snp_start (0-based)
        9: snp_end   (0-based, half-open)
        10: snp_id
        11: snp_log10p

Output:
    - <prefix>_motif_snp_summary.tsv       (one row per SNP–motif pair, BP/TP presence & scores)
    - <prefix>_motif_snp_TP_only.tsv       (subset where motif only in TP)
    - <prefix>_motif_snp_TP_only.bed       (BED file of TP-only motifs)
    - <prefix>_motif_snp_BP_only.tsv       (subset where motif only in BP)
    - <prefix>_motif_snp_BP_only.bed       (BED file of BP-only motifs)
    - <prefix>_motif_snp_both.tsv          (subset where motif in both BP and TP)
    - <prefix>_motif_snp_both.bed          (BED file of shared motifs)
    - Print summary stats to stdout.

Usage:
    python analyze_motif_snp_overlap.py \
        --input SPTBN4_motif_SNP_BP_TP_overlap.bed \
        --output_prefix SPTBN4
"""

import argparse
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Summarize BP vs TP motif differences at SNPs.")
    p.add_argument(
        "--input", "-i",
        required=True,
        help="Input BED/TSV from bedtools intersect (motif × SNP × BP/TP)."
    )
    p.add_argument(
        "--output_prefix", "-o",
        required=True,
        help="Prefix for output TSV files (e.g., 'SPTBN4')."
    )
    return p.parse_args()


def load_data(path: str) -> pd.DataFrame:
    col_names = [
        "motif_chr", "motif_start", "motif_end",
        "motif_id", "fimo_score", "hap", "fimo_pval",
        "snp_chr", "snp_start", "snp_end",
        "snp_id", "snp_log10p",
    ]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=col_names)
    return df


def summarize_motifs(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each (snp_id, motif_id), check presence in BP vs TP haplotypes.
    Keep separate coordinates for each haplotype to properly identify shared motifs
    even when they occur at slightly different positions.
    """
    # Normalize hap labels just in case
    df["hap"] = df["hap"].astype(str).str.upper()

    # Aggregate: best score per SNP/motif/hap
    # Also keep motif coordinates (from first occurrence with max score)
    grouped = (
        df.groupby(["snp_id", "motif_id", "hap"])
        .agg(
            best_score=("fimo_score", "max"),
            motif_chr=("motif_chr", "first"),
            motif_start=("motif_start", "first"),
            motif_end=("motif_end", "first"),
            snp_chr=("snp_chr", "first"),
            snp_start=("snp_start", "first"),
            snp_end=("snp_end", "first"),
            snp_log10p=("snp_log10p", "first"),
        )
        .reset_index()
    )

    # Separate BP and TP data
    bp_data = grouped[grouped["hap"] == "BP"].copy()
    tp_data = grouped[grouped["hap"] == "TP"].copy()

    # Rename columns for BP
    bp_data = bp_data.rename(columns={
        "best_score": "BP_score",
        "motif_start": "BP_motif_start",
        "motif_end": "BP_motif_end",
    })
    bp_data = bp_data.drop(columns=["hap"])

    # Rename columns for TP
    tp_data = tp_data.rename(columns={
        "best_score": "TP_score",
        "motif_start": "TP_motif_start",
        "motif_end": "TP_motif_end",
    })
    tp_data = tp_data.drop(columns=["hap"])

    # Merge BP and TP data on snp_id, motif_id (NOT on motif coordinates)
    # This ensures the same motif at different positions is treated as shared
    summary = pd.merge(
        bp_data,
        tp_data,
        on=["snp_id", "motif_id", "motif_chr", "snp_chr", "snp_start", "snp_end", "snp_log10p"],
        how="outer",
    )

    # Add presence flags
    summary["BP_present"] = summary["BP_score"].notna()
    summary["TP_present"] = summary["TP_score"].notna()

    # Add score difference (TP - BP) where both present
    summary["score_delta_TP_minus_BP"] = summary.apply(
        lambda row: (
            row["TP_score"] - row["BP_score"]
            if row["BP_present"] and row["TP_present"]
            else pd.NA
        ),
        axis=1,
    )

    # Clean up types
    summary["snp_log10p"] = pd.to_numeric(summary["snp_log10p"], errors="coerce")

    # Sort: strongest SNPs first, then motif
    summary = summary.sort_values(
        ["snp_log10p", "snp_id", "motif_id"],
        ascending=[False, True, True]
    )

    # Reorder columns for readability
    cols = [
        "snp_id", "snp_chr", "snp_start", "snp_end", "snp_log10p",
        "motif_id", "motif_chr",
        "BP_present", "BP_score", "BP_motif_start", "BP_motif_end",
        "TP_present", "TP_score", "TP_motif_start", "TP_motif_end",
        "score_delta_TP_minus_BP",
    ]
    summary = summary[cols]

    return summary


def write_bed_file(df: pd.DataFrame, output_path: str, score_column: str, haplotype: str = None):
    """
    Write a BED file from the motif summary dataframe.

    Args:
        df: DataFrame with motif data (must have motif_chr, BP_motif_start/end, TP_motif_start/end)
        output_path: Path to write BED file
        score_column: Column to use for BED score field (e.g., 'BP_score', 'TP_score', 'snp_log10p')
        haplotype: 'BP', 'TP', or 'BOTH' - determines which coordinates to use
    """
    if df.empty:
        # Write empty file
        with open(output_path, 'w') as f:
            pass
        return

    bed = pd.DataFrame()
    bed["chrom"] = df["motif_chr"]

    # Use appropriate coordinates based on haplotype
    if haplotype == "BP":
        bed["start"] = df["BP_motif_start"].astype(int)
        bed["end"] = df["BP_motif_end"].astype(int)
    elif haplotype == "TP":
        bed["start"] = df["TP_motif_start"].astype(int)
        bed["end"] = df["TP_motif_end"].astype(int)
    elif haplotype == "BOTH":
        # For shared motifs, use BP coordinates as reference
        # (could also use TP or the union of both regions)
        bed["start"] = df["BP_motif_start"].astype(int)
        bed["end"] = df["BP_motif_end"].astype(int)
    else:
        # Default: use BP if present, otherwise TP
        bed["start"] = df["BP_motif_start"].fillna(df["TP_motif_start"]).astype(int)
        bed["end"] = df["BP_motif_end"].fillna(df["TP_motif_end"]).astype(int)

    # Name: motif_id + SNP info
    if haplotype and haplotype != "BOTH":
        bed["name"] = df["motif_id"] + "_" + df["snp_id"] + "_" + haplotype
    else:
        bed["name"] = df["motif_id"] + "_" + df["snp_id"]

    # Score: use specified column, convert to int (BED spec), handle NAs
    if score_column in df.columns:
        # For FIMO scores (0-1000 range) or log10p values
        bed["score"] = df[score_column].fillna(0).apply(lambda x: int(min(1000, max(0, x * 100))))
    else:
        bed["score"] = 0

    bed["strand"] = "."

    # Write BED file (no header, tab-separated)
    bed.to_csv(output_path, sep="\t", header=False, index=False)


def main():
    args = parse_args()
    df = load_data(args.input)

    if df.empty:
        print(f"[ERROR] Input file {args.input} is empty or could not be parsed.", file=sys.stderr)
        sys.exit(1)

    summary = summarize_motifs(df)

    # Write main summary
    out_all = f"{args.output_prefix}_motif_snp_summary.tsv"
    summary.to_csv(out_all, sep="\t", index=False)

    # TP-only: present in TP, absent in BP
    tp_only = summary[(summary["TP_present"]) & (~summary["BP_present"])]
    out_tp_only = f"{args.output_prefix}_motif_snp_TP_only.tsv"
    tp_only.to_csv(out_tp_only, sep="\t", index=False)
    out_tp_only_bed = f"{args.output_prefix}_motif_snp_TP_only.bed"
    write_bed_file(tp_only, out_tp_only_bed, "TP_score", "TP")

    # BP-only: present in BP, absent in TP
    bp_only = summary[(summary["BP_present"]) & (~summary["TP_present"])]
    out_bp_only = f"{args.output_prefix}_motif_snp_BP_only.tsv"
    bp_only.to_csv(out_bp_only, sep="\t", index=False)
    out_bp_only_bed = f"{args.output_prefix}_motif_snp_BP_only.bed"
    write_bed_file(bp_only, out_bp_only_bed, "BP_score", "BP")

    # Both: motif present in both BP and TP
    both = summary[(summary["BP_present"]) & (summary["TP_present"])]
    out_both = f"{args.output_prefix}_motif_snp_both.tsv"
    both.to_csv(out_both, sep="\t", index=False)
    out_both_bed = f"{args.output_prefix}_motif_snp_both.bed"
    write_bed_file(both, out_both_bed, "snp_log10p", "BOTH")

    # Print some quick stats
    n_total = len(summary)
    n_tp_only = len(tp_only)
    n_bp_only = len(bp_only)
    n_both = len(both)

    print(f"[INFO] Wrote full summary to: {out_all}")
    print(f"[INFO] Wrote TP-only motifs to: {out_tp_only}")
    print(f"[INFO] Wrote TP-only motifs BED to: {out_tp_only_bed}")
    print(f"[INFO] Wrote BP-only motifs to: {out_bp_only}")
    print(f"[INFO] Wrote BP-only motifs BED to: {out_bp_only_bed}")
    print(f"[INFO] Wrote BP+TP motifs to: {out_both}")
    print(f"[INFO] Wrote BP+TP motifs BED to: {out_both_bed}")
    print()
    print(f"[SUMMARY] Total SNP–motif pairs: {n_total}")
    print(f"[SUMMARY] TP-only: {n_tp_only}")
    print(f"[SUMMARY] BP-only: {n_bp_only}")
    print(f"[SUMMARY] Present in both BP and TP: {n_both}")

    # Optional: how many of the "both" have non-zero score deltas?
    nonzero_delta = both["score_delta_TP_minus_BP"].dropna()
    n_nonzero = (nonzero_delta != 0).sum()
    print(f"[SUMMARY] SNP–motif pairs with motif in both BP and TP:")
    print(f"          {len(nonzero_delta)} pairs have scores in both")
    print(f"          {n_nonzero} pairs with non-zero (TP - BP) score")


if __name__ == "__main__":
    main()
