#!/usr/bin/env python3
"""
Compare GWAS peak TFBS statistics to empirical null from background windows.

For each peak locus (sptbn4b, cntnap5, itih3), loads the per-window
*_motif_snp_summary.tsv files from background windows and computes the same
BP_only_fraction statistic. Compares the observed peak value to the null
distribution and reports per-locus empirical p-values and z-scores.

Also does a per-TF-family comparison: for BCL/CEBP and other highlighted
families, tests whether peak BP-only counts exceed the background.

Usage:
  python compare_peak_to_background.py \\
    --manifest output_data/15_TFBS_Scan/background/background_windows_manifest.tsv \\
    --peak_dir output_data/15_TFBS_Scan \\
    --background_dir output_data/15_TFBS_Scan/background \\
    --output_dir output_data/15_TFBS_Scan
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


LOCI = ["sptbn4b", "cntnap5", "itih3"]

MIN_BG_WINDOWS_FOR_TEST = 10  # minimum background windows containing a TF to include it in FDR


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. Returns q-values same length as pvals."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    qvals = np.empty(n)
    qvals[order] = pvals[order] * n / (np.arange(n) + 1)
    # enforce monotonicity (cumulative min from the right)
    qvals[order] = np.minimum.accumulate(qvals[order][::-1])[::-1]
    return np.clip(qvals, 0, 1)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare GWAS peak TFBS statistics to background window null."
    )
    p.add_argument(
        "--manifest", required=True,
        help="background_windows_manifest.tsv from sample_background_windows.py",
    )
    p.add_argument(
        "--peak_dir", required=True,
        help="Directory containing per-locus peak results (output_data/15_TFBS_Scan/)",
    )
    p.add_argument(
        "--background_dir", required=True,
        help="Directory containing per-window background results (…/background/)",
    )
    p.add_argument(
        "--output_dir", default=None,
        help="Directory for output files (default: same as --background_dir)",
    )
    p.add_argument(
        "--min_snps", type=int, default=3,
        help="Minimum SNPs in a background window to include it (default: 3)",
    )
    return p.parse_args()


def load_summary(path: Path) -> pd.DataFrame:
    """Load a *_motif_snp_summary.tsv, coercing BP_present / TP_present to bool."""
    df = pd.read_csv(path, sep="\t")
    for col in ["BP_present", "TP_present"]:
        if col in df.columns and df[col].dtype != bool:
            df[col] = df[col].astype(str).str.lower().map(
                {"true": True, "false": False, "1": True, "0": False}
            )
    return df


def window_stats(df: pd.DataFrame) -> dict:
    """Compute per-window statistics from a motif_snp_summary DataFrame."""
    n_bp = int(((df["BP_present"]) & (~df["TP_present"])).sum())
    n_tp = int(((~df["BP_present"]) & (df["TP_present"])).sum())
    n_both = int(((df["BP_present"]) & (df["TP_present"])).sum())
    # Count by position in case snp_id is a placeholder (background windows use chr_pos IDs)
    if "snp_start" in df.columns:
        n_snps = int(df[["snp_chr", "snp_start"]].drop_duplicates().shape[0])
    elif "snp_id" in df.columns:
        n_snps = int(df["snp_id"].nunique())
    else:
        n_snps = 0
    n_dir = n_bp + n_tp
    bp_frac = n_bp / n_dir if n_dir > 0 else np.nan
    bp_per_snp = n_bp / n_snps if n_snps > 0 else np.nan
    return {
        "n_BP_only": n_bp,
        "n_TP_only": n_tp,
        "n_both": n_both,
        "n_snps": n_snps,
        "n_directional": n_dir,
        "BP_only_fraction": bp_frac,
        "BP_only_per_snp": bp_per_snp,
    }


def extract_tf_name(motif_id: str) -> str:
    """Extract TF name from JASPAR motif ID (MA####.#_TFNAME → TFNAME)."""
    if "_" in motif_id:
        return motif_id.split("_", 1)[1]
    return motif_id


def family_stats_from_df(df: pd.DataFrame) -> pd.DataFrame:
    """Per-TF-family BP_only and TP_only counts from a summary DataFrame."""
    df = df.copy()
    df["tf_name"] = df["motif_id"].apply(extract_tf_name)
    bp_mask = (df["BP_present"]) & (~df["TP_present"])
    tp_mask = (~df["BP_present"]) & (df["TP_present"])
    result = (
        df.groupby("tf_name")
        .apply(
            lambda g: pd.Series({
                "n_BP_only": bp_mask.loc[g.index].sum(),
                "n_TP_only": tp_mask.loc[g.index].sum(),
            }),
            include_groups=False,
        )
        .reset_index()
    )
    return result


def plot_null_distributions(locus_results: dict, output_path: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available; skipping plot.", file=sys.stderr)
        return

    n_loci = len(locus_results)
    fig, axes = plt.subplots(1, n_loci, figsize=(5 * n_loci, 4), sharey=False)
    if n_loci == 1:
        axes = [axes]

    for ax, (locus, res) in zip(axes, locus_results.items()):
        null = np.array(res["null_fractions"])
        observed = res["observed_bp_frac"]
        emp_p = res["empirical_p"]

        null_valid = null[~np.isnan(null)]
        ax.hist(null_valid, bins=30, color="steelblue", alpha=0.7, edgecolor="white", label="Background")
        if not np.isnan(observed):
            ax.axvline(observed, color="crimson", linewidth=2,
                       label=f"Peak ({observed:.3f})")
        ax.axvline(0.5, color="black", linewidth=1, linestyle="--", label="50%")
        ax.set_title(f"{locus}\np={emp_p:.3f}" if not np.isnan(emp_p) else locus, fontsize=11)
        ax.set_xlabel("BP-only fraction", fontsize=10)
        ax.set_ylabel("Background windows", fontsize=10)
        ax.legend(fontsize=8)

    fig.suptitle("GWAS peak vs. background null: BP-only motif fraction", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved null distribution plot to: {output_path}")


def main():
    args = parse_args()
    peak_dir = Path(args.peak_dir)
    bg_dir = Path(args.background_dir)
    out_dir = Path(args.output_dir) if args.output_dir else bg_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t")

    locus_results = {}
    all_family_rows = []
    locus_level_rows = []
    all_peak_vs_bg_frames = []

    print("=" * 60)
    print("TFBS Peak vs. Background Comparison")
    print("=" * 60)

    for locus in LOCI:
        print(f"\n── {locus} ──")

        # Load peak statistics
        peak_tsv = peak_dir / locus / f"{locus}_motif_snp_summary.tsv"
        if not peak_tsv.exists():
            print(f"[WARN] Peak summary not found: {peak_tsv}. Skipping.", file=sys.stderr)
            continue
        peak_df = load_summary(peak_tsv)
        peak_s = window_stats(peak_df)
        print(f"  Peak  : BP_only={peak_s['n_BP_only']}, TP_only={peak_s['n_TP_only']}, "
              f"n_SNPs={peak_s['n_snps']}, BP_frac={peak_s['BP_only_fraction']:.3f}")

        # Load background window summaries for this locus
        locus_manifest = manifest[manifest["locus"] == locus]
        if locus_manifest.empty:
            print(f"[WARN] No background windows in manifest for {locus}.", file=sys.stderr)
            continue

        bg_stats = []
        n_missing = 0
        for _, row in locus_manifest.iterrows():
            wid = row["window_id"]
            bg_tsv = bg_dir / wid / f"{wid}_motif_snp_summary.tsv"
            if not bg_tsv.exists():
                n_missing += 1
                continue
            try:
                df = load_summary(bg_tsv)
                s = window_stats(df)
                s["window_id"] = wid
                if s["n_snps"] >= args.min_snps:
                    bg_stats.append(s)
            except Exception as e:
                print(f"[WARN] Could not parse {bg_tsv}: {e}", file=sys.stderr)

        if n_missing > 0:
            print(f"  [WARN] {n_missing}/{len(locus_manifest)} background result files missing "
                  f"(run run_background_tfbs.sh first)", file=sys.stderr)

        if not bg_stats:
            print(f"  [WARN] No usable background windows for {locus} (need ≥{args.min_snps} SNPs).", file=sys.stderr)
            continue

        bg_df = pd.DataFrame(bg_stats)
        null_fracs = bg_df["BP_only_fraction"].dropna().to_numpy()

        observed = peak_s["BP_only_fraction"]
        if not np.isnan(observed) and len(null_fracs) > 0:
            null_mean = null_fracs.mean()
            null_std = null_fracs.std()
            z_score = (observed - null_mean) / null_std if null_std > 0 else np.nan
            # Two-tailed: fraction of background windows at least as extreme in either direction
            deviation = abs(observed - null_mean)
            emp_p = (np.abs(null_fracs - null_mean) >= deviation).mean()
            direction = "BP-biased" if observed > null_mean else "TP-biased"
        else:
            emp_p = np.nan
            z_score = np.nan
            direction = "NA"

        print(f"  Background ({len(null_fracs)} windows): "
              f"mean={null_fracs.mean():.3f} ± {null_fracs.std():.3f}")
        print(f"  Direction         : {direction}")
        print(f"  Empirical p-value : {emp_p:.4f}  (two-tailed)")
        print(f"  Z-score           : {z_score:.3f}" if not np.isnan(z_score) else "  Z-score           : NA")

        locus_results[locus] = {
            "peak_stats": peak_s,
            "null_fractions": null_fracs.tolist(),
            "observed_bp_frac": observed,
            "empirical_p": emp_p,
            "z_score": z_score,
            "direction": direction,
            "n_bg_windows": len(null_fracs),
            "sig_tfs": [],  # populated below after per-TF analysis
        }

        locus_level_rows.append({
            "locus": locus,
            "n_snps": peak_s["n_snps"],
            "n_bp_only": peak_s["n_BP_only"],
            "n_tp_only": peak_s["n_TP_only"],
            "bp_only_fraction": round(observed, 4) if not np.isnan(observed) else np.nan,
            "n_bg_windows": len(null_fracs),
            "bg_mean_bp_only_fraction": round(null_fracs.mean(), 4),
            "bg_sd_bp_only_fraction": round(null_fracs.std(), 4),
            "z_score": round(z_score, 3) if not np.isnan(z_score) else np.nan,
            "empirical_p_twotailed": emp_p,
            "direction": direction,
            "n_sig_tfs_q0.20": 0,  # filled after per-TF section
        })

        # Per-family comparison: peak counts vs background distribution
        peak_families = family_stats_from_df(peak_df)
        peak_families["locus"] = locus
        peak_families["source"] = "peak"
        all_family_rows.append(peak_families)

        bg_family_frames = []
        for _, row in locus_manifest.iterrows():
            wid = row["window_id"]
            bg_tsv = bg_dir / wid / f"{wid}_motif_snp_summary.tsv"
            if not bg_tsv.exists():
                continue
            try:
                df = load_summary(bg_tsv)
                s = window_stats(df)
                if s["n_snps"] < args.min_snps:
                    continue
                fam = family_stats_from_df(df)
                fam["window_id"] = wid
                bg_family_frames.append(fam)
            except Exception:
                pass

        if bg_family_frames:
            bg_family_df = pd.concat(bg_family_frames, ignore_index=True)
            n_bg = len(bg_family_frames)

            # Per-TF background summary: mean, SD, and window count
            bg_agg = (
                bg_family_df.groupby("tf_name")["n_BP_only"]
                .agg(bg_mean_BP_only="mean", bg_sd_BP_only="std", n_bg_windows_with_tf="count")
                .reset_index()
            )
            bg_tp_mean = (
                bg_family_df.groupby("tf_name")["n_TP_only"]
                .mean()
                .rename("bg_mean_TP_only")
                .reset_index()
            )
            bg_agg = bg_agg.merge(bg_tp_mean, on="tf_name", how="left")

            peak_vs_bg = peak_families.merge(bg_agg, on="tf_name", how="left")
            peak_vs_bg["locus"] = locus
            peak_vs_bg["n_bg_windows"] = n_bg

            # Two-tailed empirical p-value per TF + direction
            emp_ps = []
            directions = []
            for _, frow in peak_vs_bg.iterrows():
                tf = frow["tf_name"]
                peak_n = frow["n_BP_only"]
                bg_counts = bg_family_df[bg_family_df["tf_name"] == tf]["n_BP_only"]
                if len(bg_counts) > 0:
                    null_mean = bg_counts.mean()
                    deviation = abs(peak_n - null_mean)
                    ep = (np.abs(bg_counts - null_mean) >= deviation).mean()
                    d = "BP-enriched" if peak_n > null_mean else "TP-enriched"
                else:
                    ep = np.nan
                    d = "NA"
                emp_ps.append(ep)
                directions.append(d)
            peak_vs_bg["direction"] = directions
            peak_vs_bg["family_empirical_p"] = emp_ps

            # BH-FDR correction — only for TFs with enough background coverage
            n_rows = len(peak_vs_bg)
            q_values = np.full(n_rows, np.nan)
            bg_coverage = peak_vs_bg["n_bg_windows_with_tf"].fillna(0).to_numpy()
            testable_pos = np.where(bg_coverage >= MIN_BG_WINDOWS_FOR_TEST)[0]
            if len(testable_pos) > 0:
                raw_p = peak_vs_bg["family_empirical_p"].to_numpy(dtype=float)[testable_pos]
                raw_p = np.where(np.isnan(raw_p), 1.0, raw_p)
                q_values[testable_pos] = _bh_fdr(raw_p)
            peak_vs_bg["q_value"] = q_values

            fam_out = out_dir / f"tfbs_family_background_{locus}.csv"
            cols = ["locus", "tf_name", "n_BP_only", "n_TP_only",
                    "bg_mean_BP_only", "bg_sd_BP_only", "bg_mean_TP_only",
                    "n_bg_windows", "n_bg_windows_with_tf",
                    "direction", "family_empirical_p", "q_value"]
            peak_vs_bg[cols].sort_values(
                ["q_value", "family_empirical_p"], na_position="last"
            ).to_csv(fam_out, index=False)
            print(f"  Family comparison written to: {fam_out}")

            # Report significant TFs to stdout
            sig = peak_vs_bg[peak_vs_bg["q_value"] < 0.20].sort_values("q_value")
            sig_tf_names = sig["tf_name"].tolist()
            locus_results[locus]["sig_tfs"] = sig_tf_names
            locus_level_rows[-1]["n_sig_tfs_q0.20"] = len(sig_tf_names)
            all_peak_vs_bg_frames.append(peak_vs_bg[cols])
            if len(sig) > 0:
                print(f"  Significant TFs (q<0.20): {len(sig)}")
                for _, r in sig.head(10).iterrows():
                    print(f"    {r['tf_name']:30s}  {r['direction']:12s}  "
                          f"peak={int(r['n_BP_only'])}  bg_mean={r['bg_mean_BP_only']:.2f}  "
                          f"p={r['family_empirical_p']:.3f}  q={r['q_value']:.3f}")
            else:
                top = peak_vs_bg.dropna(subset=["family_empirical_p"]).sort_values("family_empirical_p").head(10)
                print(f"  No TFs with q<0.20. Top 10 by raw p-value:")
                for _, r in top.iterrows():
                    print(f"    {r['tf_name']:30s}  {r['direction']:12s}  "
                          f"peak={int(r['n_BP_only'])}  bg_mean={r['bg_mean_BP_only']:.2f}  "
                          f"p={r['family_empirical_p']:.3f}")

    # ── Summary text file ──────────────────────────────────────────────────────
    summary_path = out_dir / "tfbs_background_comparison.txt"
    with open(summary_path, "w") as fh:
        fh.write("TFBS Peak vs. Background Window Comparison\n")
        fh.write("=" * 50 + "\n\n")
        for locus, res in locus_results.items():
            ps = res["peak_stats"]
            fh.write(f"Locus: {locus}\n")
            fh.write(f"  Peak BP_only: {ps['n_BP_only']}, TP_only: {ps['n_TP_only']}, "
                     f"n_SNPs: {ps['n_snps']}\n")
            fh.write(f"  Observed BP_only_fraction: {res['observed_bp_frac']:.4f}\n")
            fh.write(f"  Background windows used: {res['n_bg_windows']}\n")
            fh.write(f"  Background mean ± SD: "
                     f"{np.mean(res['null_fractions']):.4f} ± {np.std(res['null_fractions']):.4f}\n")
            fh.write(f"  Direction: {res['direction']}\n")
            fh.write(f"  Empirical p-value (two-tailed): {res['empirical_p']:.4f}\n")
            fh.write(f"  Z-score: {res['z_score']:.3f}\n")
            sig = res.get("sig_tfs", [])
            if sig:
                fh.write(f"  Significant TFs (q<0.20, n={len(sig)}): {', '.join(sig)}\n")
            else:
                fh.write(f"  Significant TFs (q<0.20): none\n")
            fh.write("\n")
    print(f"\n[INFO] Summary written to: {summary_path}")

    # ── Supplemental TSVs ──────────────────────────────────────────────────────
    if locus_level_rows:
        locus_tsv = out_dir / "tfbs_locus_level_summary.tsv"
        pd.DataFrame(locus_level_rows).to_csv(locus_tsv, sep="\t", index=False)
        print(f"[INFO] Locus-level summary written to: {locus_tsv}")

    if all_peak_vs_bg_frames:
        tf_tsv = out_dir / "tfbs_per_tf_enrichment.tsv"
        tf_df = pd.concat(all_peak_vs_bg_frames, ignore_index=True)
        tf_df["testable"] = tf_df["n_bg_windows_with_tf"].fillna(0) >= MIN_BG_WINDOWS_FOR_TEST
        tf_df["significant_q0.20"] = tf_df["q_value"] < 0.20
        tf_df.sort_values(["locus", "q_value", "family_empirical_p"], na_position="last").to_csv(
            tf_tsv, sep="\t", index=False
        )
        print(f"[INFO] Per-TF enrichment table written to: {tf_tsv}")

    # ── Annotate motif_snp_summary.tsv files with statistical columns ─────────
    if all_peak_vs_bg_frames and locus_results:
        tf_enrich = pd.concat(all_peak_vs_bg_frames, ignore_index=True)
        tf_enrich["tf_testable"] = tf_enrich["n_bg_windows_with_tf"].fillna(0) >= MIN_BG_WINDOWS_FOR_TEST
        tf_enrich["tf_significant_q0.20"] = tf_enrich["q_value"] < 0.20

        for locus, res in locus_results.items():
            peak_tsv = peak_dir / locus / f"{locus}_motif_snp_summary.tsv"
            if not peak_tsv.exists():
                continue
            df = load_summary(peak_tsv)
            df["tf_name"] = df["motif_id"].apply(extract_tf_name)

            # Level 1: locus-wide empirical test (same value on every row)
            df["locus_bp_only_fraction"] = round(res["observed_bp_frac"], 4)
            df["locus_bg_mean_bp_fraction"] = round(float(np.mean(res["null_fractions"])), 4)
            df["locus_empirical_p"] = res["empirical_p"]
            df["locus_direction"] = res["direction"]

            # Level 2: per-TF enrichment stats (one lookup per unique TF name)
            tf_locus = (
                tf_enrich[tf_enrich["locus"] == locus][
                    ["tf_name", "bg_mean_BP_only", "direction",
                     "family_empirical_p", "q_value",
                     "tf_testable", "tf_significant_q0.20"]
                ]
                .drop_duplicates("tf_name")
                .rename(columns={
                    "bg_mean_BP_only": "tf_bg_mean_BP_only",
                    "direction": "tf_direction",
                    "family_empirical_p": "tf_family_empirical_p",
                    "q_value": "tf_q_value",
                })
            )
            df = df.merge(tf_locus, on="tf_name", how="left")
            df.to_csv(peak_tsv, sep="\t", index=False)
            print(f"[INFO] Annotated: {peak_tsv}")

    # ── Plot ───────────────────────────────────────────────────────────────────
    if locus_results:
        plot_path = out_dir / "tfbs_background_null_distributions.png"
        plot_null_distributions(locus_results, plot_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
