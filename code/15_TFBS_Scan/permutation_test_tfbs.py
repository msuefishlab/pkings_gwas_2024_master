#!/usr/bin/env python3
"""
Permutation test and binomial family test for TFBS directionality.

Addresses reviewer concern that BP-biased motif gains at GWAS peaks may not be
statistically significant compared to a null background.

Two tests:
  1. Allele-label permutation test: is the overall BP_only_fraction significantly
     greater than expected if BP/TP labels were random?
  2. Per-TF-family binomial test: within each TF family, is the number of BP-only
     hits significantly greater than 50% of (BP-only + TP-only)?

Usage:
  python permutation_test_tfbs.py \\
    -i SPTBN4B:output_data/15_TFBS_Scan/sptbn4b/sptbn4b_motif_snp_summary.tsv \\
    -i CNTNAP5:output_data/15_TFBS_Scan/cntnap5/cntnap5_motif_snp_summary.tsv \\
    -i ITIH3:output_data/15_TFBS_Scan/itih3/itih3_motif_snp_summary.tsv \\
    --output_dir output_data/15_TFBS_Scan \\
    --n_permutations 10000 \\
    --seed 42
"""

import argparse
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from scipy.stats import binomtest


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Permutation and binomial tests for BP-biased TFBS directionality."
    )
    p.add_argument(
        "-i", "--input",
        action="append",
        required=True,
        metavar="GENE:PATH",
        help="Input TSV in the form GENE_LABEL:/path/to/file.tsv. Can be given multiple times.",
    )
    p.add_argument(
        "--output_dir", "-o",
        default=".",
        help="Directory for output files (default: current directory).",
    )
    p.add_argument(
        "--n_permutations", "-n",
        type=int,
        default=10000,
        help="Number of permutations for Test 1 (default: 10000).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    p.add_argument(
        "--min_family_n",
        type=int,
        default=3,
        help="Minimum (n_BP_only + n_TP_only) per TF family to include in Test 2 (default: 3).",
    )
    return p.parse_args()


def parse_gene_inputs(raw_inputs) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for item in raw_inputs:
        if ":" not in item:
            raise ValueError(f"Input '{item}' must be in form GENE:PATH")
        gene, path = item.split(":", 1)
        gene, path = gene.strip(), path.strip()
        if gene in mapping:
            raise ValueError(f"Duplicate gene label '{gene}'")
        mapping[gene] = path
    return mapping


def load_data(gene_to_path: Dict[str, str]) -> pd.DataFrame:
    frames = []
    for gene, path in gene_to_path.items():
        df = pd.read_csv(path, sep="\t")
        df["locus"] = gene
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)

    for col in ["BP_present", "TP_present"]:
        if combined[col].dtype != bool:
            combined[col] = combined[col].astype(str).str.lower().map(
                {"true": True, "false": False, "1": True, "0": False}
            )
    return combined


def bp_only_fraction(df: pd.DataFrame) -> float:
    """Fraction of directional motif-SNP pairs that are BP-only."""
    n_bp = ((df["BP_present"]) & (~df["TP_present"])).sum()
    n_tp = ((~df["BP_present"]) & (df["TP_present"])).sum()
    total = n_bp + n_tp
    if total == 0:
        return np.nan
    return n_bp / total


def run_permutation_test(df: pd.DataFrame, n_permutations: int, seed: int):
    """
    Permute BP/TP labels at the SNP level and rebuild the BP_only_fraction.

    For each permutation, each unique SNP is independently assigned a coin flip:
    heads → keep labels, tails → swap BP_present and TP_present for all motifs
    at that SNP. This preserves the multi-motif correlation structure per SNP.
    """
    rng = np.random.default_rng(seed)
    observed = bp_only_fraction(df)

    snp_ids = df["snp_id"].unique()
    n_snps = len(snp_ids)
    snp_index = {snp: i for i, snp in enumerate(snp_ids)}
    snp_col = df["snp_id"].map(snp_index).to_numpy()

    bp_arr = df["BP_present"].to_numpy(dtype=bool)
    tp_arr = df["TP_present"].to_numpy(dtype=bool)

    null_fractions = np.empty(n_permutations)

    for i in range(n_permutations):
        swap = rng.integers(0, 2, size=n_snps, dtype=bool)  # True = swap this SNP
        flip_mask = swap[snp_col]

        perm_bp = np.where(flip_mask, tp_arr, bp_arr)
        perm_tp = np.where(flip_mask, bp_arr, tp_arr)

        n_bp = ((perm_bp) & (~perm_tp)).sum()
        n_tp = ((~perm_bp) & (perm_tp)).sum()
        total = n_bp + n_tp
        null_fractions[i] = n_bp / total if total > 0 else np.nan

    empirical_p = (null_fractions >= observed).mean()
    return observed, null_fractions, empirical_p


def extract_tf_name(motif_id: str) -> str:
    """
    Extract the TF name from a JASPAR motif ID.
    JASPAR format: MA####.#_TFNAME  e.g. MA0139.1_CTCF → CTCF
    Falls back to the full motif_id if the pattern is not matched.
    """
    if "_" in motif_id:
        return motif_id.split("_", 1)[1]
    return motif_id


def run_family_binomial_tests(df: pd.DataFrame, min_n: int) -> pd.DataFrame:
    """
    For each TF family, test if n_BP_only is significantly > 50% of directional pairs.
    Uses a one-sided binomial test (alternative='greater').
    Applies Benjamini-Hochberg FDR correction.
    """
    df = df.copy()
    df["tf_name"] = df["motif_id"].apply(extract_tf_name)

    bp_only_mask = (df["BP_present"]) & (~df["TP_present"])
    tp_only_mask = (~df["BP_present"]) & (df["TP_present"])

    family_stats = (
        df.groupby("tf_name")
        .apply(
            lambda g: pd.Series({
                "n_BP_only": bp_only_mask.loc[g.index].sum(),
                "n_TP_only": tp_only_mask.loc[g.index].sum(),
                "n_loci": g["locus"].nunique(),
            }),
            include_groups=False,
        )
        .reset_index()
    )

    family_stats["n_directional"] = family_stats["n_BP_only"] + family_stats["n_TP_only"]
    family_stats = family_stats[family_stats["n_directional"] >= min_n].copy()

    if family_stats.empty:
        print(f"[WARN] No TF families with n_directional >= {min_n}.", file=sys.stderr)
        family_stats["p_value"] = np.nan
        family_stats["q_value"] = np.nan
        return family_stats

    family_stats["p_value"] = family_stats.apply(
        lambda r: binomtest(
            int(r["n_BP_only"]),
            int(r["n_directional"]),
            p=0.5,
            alternative="greater",
        ).pvalue,
        axis=1,
    )

    family_stats["q_value"] = _bh_fdr(family_stats["p_value"].to_numpy())
    family_stats = family_stats.sort_values("p_value")

    return family_stats


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (no external dependency required)."""
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    qvals = np.minimum(1.0, pvals * n / ranks)
    # Enforce monotonicity from the right
    qvals_ordered = qvals[order]
    for i in range(n - 2, -1, -1):
        qvals_ordered[i] = min(qvals_ordered[i], qvals_ordered[i + 1])
    qvals[order] = qvals_ordered
    return qvals


def plot_null_distribution(null_fractions, observed, output_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available; skipping null distribution plot.", file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(null_fractions, bins=50, color="steelblue", alpha=0.75, edgecolor="white",
            label="Null (permuted)")
    ax.axvline(observed, color="crimson", linewidth=2, label=f"Observed ({observed:.3f})")
    ax.axvline(0.5, color="black", linewidth=1, linestyle="--", label="50% (no bias)")
    ax.set_xlabel("BP-only fraction", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Permutation null distribution: BP-only fraction", fontsize=13)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"[INFO] Saved null distribution plot to: {output_path}")


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_to_path = parse_gene_inputs(args.input)
    df = load_data(gene_to_path)

    n_total = len(df)
    n_bp_only = ((df["BP_present"]) & (~df["TP_present"])).sum()
    n_tp_only = ((~df["BP_present"]) & (df["TP_present"])).sum()
    n_both = ((df["BP_present"]) & (df["TP_present"])).sum()

    print("=" * 60)
    print("TFBS Directionality Statistical Tests")
    print("=" * 60)
    print(f"\nInput loci: {', '.join(gene_to_path)}")
    print(f"Total SNP-motif pairs: {n_total}")
    print(f"  BP-only : {n_bp_only}")
    print(f"  TP-only : {n_tp_only}")
    print(f"  Both    : {n_both}")

    # ── Test 1: Permutation ────────────────────────────────────────────────────
    print(f"\n── Test 1: Allele-label permutation (n={args.n_permutations:,}) ──")
    observed, null_fractions, empirical_p = run_permutation_test(
        df, args.n_permutations, args.seed
    )
    print(f"Observed BP-only fraction : {observed:.4f}")
    print(f"Null mean ± SD            : {null_fractions.mean():.4f} ± {null_fractions.std():.4f}")
    print(f"Empirical p-value         : {empirical_p:.4f}")

    plot_path = out_dir / "permutation_null_distribution.png"
    plot_null_distribution(null_fractions, observed, plot_path)

    perm_results_path = out_dir / "tfbs_permutation_results.txt"
    with open(perm_results_path, "w") as fh:
        fh.write("TFBS Allele-label Permutation Test\n")
        fh.write("=" * 40 + "\n")
        fh.write(f"Loci: {', '.join(gene_to_path)}\n")
        fh.write(f"n_permutations: {args.n_permutations}\n")
        fh.write(f"random_seed: {args.seed}\n\n")
        fh.write(f"n_total_pairs: {n_total}\n")
        fh.write(f"n_BP_only: {n_bp_only}\n")
        fh.write(f"n_TP_only: {n_tp_only}\n")
        fh.write(f"n_both: {n_both}\n\n")
        fh.write(f"observed_BP_only_fraction: {observed:.6f}\n")
        fh.write(f"null_mean: {null_fractions.mean():.6f}\n")
        fh.write(f"null_sd: {null_fractions.std():.6f}\n")
        fh.write(f"empirical_p_value: {empirical_p:.6f}\n")
    print(f"[INFO] Wrote permutation results to: {perm_results_path}")

    # ── Test 2: Per-family binomial tests ─────────────────────────────────────
    print(f"\n── Test 2: Per-TF-family binomial tests (min_n={args.min_family_n}) ──")
    family_results = run_family_binomial_tests(df, args.min_family_n)

    if not family_results.empty:
        n_sig = (family_results["q_value"] < 0.05).sum()
        print(f"TF families tested: {len(family_results)}")
        print(f"Significant after FDR (q < 0.05): {n_sig}")
        print("\nTop families by p-value:")
        cols = ["tf_name", "n_BP_only", "n_TP_only", "n_directional", "n_loci", "p_value", "q_value"]
        print(family_results[cols].head(20).to_string(index=False))

        family_out = out_dir / "tfbs_family_binomial_tests.csv"
        family_results[cols].to_csv(family_out, index=False)
        print(f"\n[INFO] Wrote family binomial results to: {family_out}")
    else:
        print("No TF families passed the minimum count filter.")

    print("\nDone.")


if __name__ == "__main__":
    main()
