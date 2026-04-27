#!/usr/bin/env python3
"""
Summarize TFBS differences between TP and BP haplotypes across multiple loci.

Expected input columns (as in your hi_signf_*_new.csv files):
  snp_id, snp_chr, snp_start, snp_end, snp_log10p,
  motif_id, motif_chr,
  BP_present, BP_score, BP_motif_start, BP_motif_end,
  TP_present, TP_score, TP_motif_start, TP_motif_end,
  score_delta_TP_minus_BP

Usage example:

  python tfbs_summary.py \
    -i SPTBN4:hi_signf_sptbn4_new.csv \
    -i ITIH3:hi_signf_itih3_new.csv \
    -i CNTNAP5:hi_signf_cntnap5_new.csv

Outputs:
  - Prints per-gene summaries to stdout.
  - Writes per-gene motif summaries:
      motif_summary_<GENE>.csv
  - Writes cross-gene motif matrix:
      motif_cross_gene_matrix.csv
"""

import argparse
import sys
from typing import Dict, Tuple
import pandas as pd
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize TFBS changes between TP and BP haplotypes."
    )
    parser.add_argument(
        "-i",
        "--input",
        action="append",
        required=True,
        metavar="GENE:PATH",
        help=(
            "Input CSV in the form GENE_LABEL:/path/to/file.csv. "
            "Can be given multiple times."
        ),
    )
    return parser.parse_args()


def parse_gene_inputs(raw_inputs) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for item in raw_inputs:
        if ":" not in item:
            raise ValueError(f"Input '{item}' must be in form GENE:PATH")
        gene, path = item.split(":", 1)
        gene = gene.strip()
        path = path.strip()
        if not gene:
            raise ValueError(f"Empty gene label in '{item}'")
        if gene in mapping:
            raise ValueError(f"Duplicate gene label '{gene}'")
        mapping[gene] = path
    return mapping


def classify_motif_row(bp_present: bool, tp_present: bool) -> str:
    """
    Classify a single row (motif–SNP entry) by presence combination.
    """
    if bp_present and tp_present:
        return "BP_and_TP"
    elif bp_present and not tp_present:
        return "BP_only"
    elif not bp_present and tp_present:
        return "TP_only"
    else:
        return "neither"


def classify_motif_any(bp_any: bool, tp_any: bool) -> str:
    """
    Classify a motif_id by whether it appears in BP and/or TP at least once.
    """
    if bp_any and tp_any:
        return "both"
    elif bp_any and not tp_any:
        return "BP_only"
    elif not bp_any and tp_any:
        return "TP_only"
    else:
        return "neither"


def summarize_gene(gene: str, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute summaries for a single gene.

    Returns:
      motif_summary (per motif_id)
      delta_stats_per_motif (per motif_id with non-NaN score_delta_TP_minus_BP)
    """

    print(f"\n========== {gene} ==========")

    # Ensure boolean dtype (in case they come in as object)
    for col in ["BP_present", "TP_present"]:
        if df[col].dtype != bool:
            # Accept True/False / 0/1 / strings
            df[col] = df[col].astype(str).str.lower().map(
                {"true": True, "false": False, "1": True, "0": False}
            )

    # ---- Row-level classification (motif–SNP entries) ----
    df["row_class"] = [
        classify_motif_row(bp, tp) for bp, tp in zip(df["BP_present"], df["TP_present"])
    ]

    row_counts = df["row_class"].value_counts().reindex(
        ["BP_only", "TP_only", "BP_and_TP", "neither"], fill_value=0
    )

    print("Row-level motif–SNP counts (based on BP_present / TP_present):")
    for cls, n in row_counts.items():
        print(f"  {cls:10s}: {n}")

    # ---- Motif-level presence (any BP / any TP) ----
    motif_group = df.groupby("motif_id", dropna=False).agg(
        BP_any=("BP_present", "any"),
        TP_any=("TP_present", "any"),
        n_rows=("motif_id", "size"),
    )
    motif_group["class"] = motif_group.apply(
        lambda r: classify_motif_any(r["BP_any"], r["TP_any"]), axis=1
    )

    motif_class_counts = motif_group["class"].value_counts().reindex(
        ["BP_only", "TP_only", "both", "neither"], fill_value=0
    )

    print("\nMotif-level counts (per motif_id):")
    for cls, n in motif_class_counts.items():
        print(f"  {cls:10s}: {n}")

    # ---- score_delta_TP_minus_BP summary ----
    if "score_delta_TP_minus_BP" in df.columns:
        delta = df["score_delta_TP_minus_BP"].dropna()
        print("\nScore delta (TP − BP) summary (per row where not NaN):")
        if len(delta) == 0:
            print("  No non-NaN score deltas for this gene.")
        else:
            n = len(delta)
            print(f"  n_rows_with_delta: {n}")
            print(f"  mean_delta      : {delta.mean():.4f}")
            print(f"  median_delta    : {delta.median():.4f}")
            print(f"  n_delta > 0     : {(delta > 0).sum()}")
            print(f"  n_delta < 0     : {(delta < 0).sum()}")
            print(f"  n_delta == 0    : {(delta == 0).sum()}")

            # motif-level delta stats (for later cross-gene matrix)
            delta_stats_per_motif = (
                df.dropna(subset=["score_delta_TP_minus_BP"])
                .groupby("motif_id")["score_delta_TP_minus_BP"]
                .agg(count="count", mean="mean", median="median")
                .reset_index()
            )
    else:
        print("\nNo 'score_delta_TP_minus_BP' column found.")
        delta_stats_per_motif = pd.DataFrame(
            columns=["motif_id", "count", "mean", "median"]
        )

    # Combine motif presence and delta stats into one per-gene motif summary
    motif_summary = motif_group.reset_index()
    if not delta_stats_per_motif.empty:
        motif_summary = motif_summary.merge(
            delta_stats_per_motif,
            on="motif_id",
            how="left",
            suffixes=("", "_delta"),
        )
    else:
        motif_summary["count"] = np.nan
        motif_summary["mean"] = np.nan
        motif_summary["median"] = np.nan

    # Save per-gene motif summary
    outfname = f"motif_summary_{gene}.csv"
    motif_summary.to_csv(outfname, index=False)
    print(f"\nWrote per-motif summary for {gene} to: {outfname}")

    return motif_summary, delta_stats_per_motif


def build_cross_gene_matrix(
    per_gene_motif_summaries: Dict[str, pd.DataFrame]
) -> pd.DataFrame:
    """
    Build a cross-gene matrix of motif classes and mean score deltas.
    """

    # Start with a union of all motif_ids
    all_motifs = sorted(
        set().union(*[set(df["motif_id"].tolist()) for df in per_gene_motif_summaries.values()])
    )
    matrix = pd.DataFrame({"motif_id": all_motifs})

    for gene, df in per_gene_motif_summaries.items():
        df_small = df[["motif_id", "class", "mean"]].copy()
        df_small.rename(
            columns={
                "class": f"class_{gene}",
                "mean": f"mean_delta_{gene}",
            },
            inplace=True,
        )
        matrix = matrix.merge(df_small, on="motif_id", how="left")

    return matrix


def main():
    args = parse_args()
    gene_to_path = parse_gene_inputs(args.input)

    per_gene_motif_summaries: Dict[str, pd.DataFrame] = {}

    for gene, path in gene_to_path.items():
        print(f"\nReading {gene} from {path}")
        df = pd.read_csv(path)
        motif_summary, _ = summarize_gene(gene, df)
        per_gene_motif_summaries[gene] = motif_summary

    # ---- Cross-gene motif matrix ----
    print("\n========== Cross-gene motif matrix ==========")
    matrix = build_cross_gene_matrix(per_gene_motif_summaries)
    outfname = "motif_cross_gene_matrix.csv"
    matrix.to_csv(outfname, index=False)
    print(f"Wrote cross-gene motif matrix to: {outfname}")

    # Optional: show motifs present in more than one gene
    cols_class = [c for c in matrix.columns if c.startswith("class_")]
    if len(cols_class) > 1:
        multi_gene = matrix.dropna(subset=cols_class, how="all").copy()
        # Motifs appearing in at least two genes (non-NaN class)
        multi_gene["n_genes"] = multi_gene[cols_class].notna().sum(axis=1)
        shared = multi_gene[multi_gene["n_genes"] >= 2]
        print(
            f"\nMotifs appearing in ≥2 genes (based on non-NaN class columns): "
            f"{len(shared)} motifs"
        )
        # Show a few as an example
        print(shared[["motif_id"] + cols_class].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
