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
    - <prefix>_motif_snp_BP_only.tsv       (subset where motif only in BP)
    - <prefix>_motif_snp_both.tsv          (subset where motif in both BP and TP)
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
    For each (snp_id, motif_id, hap), compute the maximum FIMO score.
    Then pivot to compare BP vs TP scores & presence.
    """
    # Normalize hap labels just in case
    df["hap"] = df["hap"].astype(str).str.upper()

    # Aggregate: best score per SNP/motif/hap
    grouped = (
        df.groupby(["snp_id", "motif_id", "hap"])
        .agg(
            best_score=("fimo_score", "max"),
            snp_chr=("snp_chr", "first"),
            snp_start=("snp_start", "first"),
            snp_end=("snp_end", "first"),
            snp_log10p=("snp_log10p", "first"),
        )
        .reset_index()
    )

    # Pivot into BP and TP columns
    pivot = grouped.pivot_table(
        index=["snp_id", "motif_id", "snp_chr", "snp_start", "snp_end", "snp_log10p"],
        columns="hap",
        values="best_score",
        aggfunc="max",
    ).reset_index()

    # Ensure columns exist (BP/TP may be missing if not present at all)
    if "BP" not in pivot.columns:
        pivot["BP"] = pd.NA
    if "TP" not in pivot.columns:
        pivot["TP"] = pd.NA

    pivot = pivot.rename(columns={"BP": "BP_score", "TP": "TP_score"})

    # Add presence flags
    pivot["BP_present"] = pivot["BP_score"].notna()
    pivot["TP_present"] = pivot["TP_score"].notna()

    # Add score difference (TP - BP) where both present
    pivot["score_delta_TP_minus_BP"] = pivot.apply(
        lambda row: (
            row["TP_score"] - row["BP_score"]
            if row["BP_present"] and row["TP_present"]
            else pd.NA
        ),
        axis=1,
    )

    # Clean up types
    pivot["snp_log10p"] = pd.to_numeric(pivot["snp_log10p"], errors="coerce")

    # Sort: strongest SNPs first, then motif
    pivot = pivot.sort_values(["snp_log10p", "snp_id", "motif_id"], ascending=[False, True, True])

    # Reorder columns for readability
    cols = [
        "snp_id", "snp_chr", "snp_start", "snp_end", "snp_log10p",
        "motif_id",
        "BP_present", "BP_score",
        "TP_present", "TP_score",
        "score_delta_TP_minus_BP",
    ]
    pivot = pivot[cols]

    return pivot


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

    # BP-only: present in BP, absent in TP
    bp_only = summary[(summary["BP_present"]) & (~summary["TP_present"])]
    out_bp_only = f"{args.output_prefix}_motif_snp_BP_only.tsv"
    bp_only.to_csv(out_bp_only, sep="\t", index=False)

    # Both: motif present in both BP and TP
    both = summary[(summary["BP_present"]) & (summary["TP_present"])]
    out_both = f"{args.output_prefix}_motif_snp_both.tsv"
    both.to_csv(out_both, sep="\t", index=False)

    # Print some quick stats
    n_total = len(summary)
    n_tp_only = len(tp_only)
    n_bp_only = len(bp_only)
    n_both = len(both)

    print(f"[INFO] Wrote full summary to: {out_all}")
    print(f"[INFO] Wrote TP-only motifs to: {out_tp_only}")
    print(f"[INFO] Wrote BP-only motifs to: {out_bp_only}")
    print(f"[INFO] Wrote BP+TP motifs to: {out_both}")
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
