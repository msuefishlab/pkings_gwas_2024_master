#!/usr/bin/env python3
"""10_checks.py
Spectrum and model checks of DESIGN.md sections 9.2 and 9.3. Convergence is checked in 08_collect_fits.py;
Poisson parameter recovery is run elsewhere.

Spectra, per pair: checks/spectra_checks<suffix>.tsv (pair, check, value, expected, pass)
  block_sum[_alt]     the 1-Mb block spectra sum to the genome-wide fitted spectrum
  total[_alt]         spectrum total + block invariant sites + invariant sites in blocks with no fitting-mask SNP
                      (06_gather.py keys block_inv by the SNP blocks, so those sites are in L_mask only; counted
                      here from sfs/chunks, reported as invariant_in_snp_free_blocks[_alt]) = L_mask
  fold_*              folding checks within and across the two variants (Deviation 1; see spectra_rows)
  from_data_dict      not run here (needs the site table)
Models, per pair x model in fits/best_fits<suffix>.tsv: checks/model_checks<suffix>.tsv
  pi_BP, pi_TP, dxy   per callable site (L_mask), theta-scaled model and data, ratio model / data
  Fst_WC              dadi Spectrum.Fst(), which is Weir & Cockerham (1984), not Hudson
  Fst_Hudson          1 - mean(pi_BP, pi_TP) / dxy (Hudson et al. 1992), from the columns above
  TajD_BP, TajD_TP    Tajima's D of each marginal spectrum, model and data, and whether the signs agree;
                      sign_check_required marks the BP side of BIK and BAM pairs (A2w, A1w), where
                      section 9.3 requires agreement
  corner_*_change     relative change of the expected count summed over a corner region when the grid
                      doubles (pts -> 2 pts): BP corner delta >= 0.95 and TP corner delta <= -0.95 for
                      unfolded fits, |delta| >= 0.95 (corner_either) for folded fits. Any change > 0.05 makes
                      the model's E "numerically unreliable" (section 9.3; 09_tail.py reads this column)
"BP" and "TP" are population 1 and 2 of the pair, so in the TP-TP controls they mean APA/BEN and DOV/DOG.
The expected spectrum is built as in 09_tail.py: func_ex at the best parameters, theta from
optimal_sfs_scaling (on the folded model for folded fits), misidentification included for unfolded fits.

--skip_numerics leaves out the doubled grid (M3 costs 5-11 s per evaluation at the default grid); the corner
columns are then NA, numerics_run False and numerically_unreliable False, so 09_tail.py applies no tag.

Usage (inside dadi_tools.sif, from the repo root):
  python3 .../10_checks.py [--suffix _x] [--chunk_suffix _y] [--pairs A1w,A2w] [--skip_numerics]
"""
import argparse
import collections
import csv
import glob
import importlib.util
import os
import pickle
import re
import time

import numpy as np
import dadi

OUT = os.path.join(os.getcwd(), os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("fitmod", os.path.join(HERE, "07_fit_models.py"))
fitmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitmod)

CORNER = 0.95
MAX_CHANGE = 0.05            # section 9.3: corner cells must move by less than 5% when the grid doubles
EPS = 1e-9                   # floating-point guard on the corner thresholds
RTOL, ATOL = 1e-8, 1e-6      # summation-order tolerance for spectra built from projected (fractional) counts
SIGN_POPS = ("BAM", "BIK")   # section 9.3: the sign of Tajima's D must match for these (BP side)
MODEL_COLS = ["pair", "model", "analysis", "converged", "L_mask",
              "pi_BP_model", "pi_BP_data", "pi_BP_ratio", "pi_TP_model", "pi_TP_data", "pi_TP_ratio",
              "dxy_model", "dxy_data", "dxy_ratio", "Fst_WC_model", "Fst_WC_data", "Fst_WC_ratio",
              "Fst_Hudson_model", "Fst_Hudson_data", "Fst_Hudson_ratio",
              "TajD_BP_model", "TajD_BP_data", "TajD_BP_sign_match",
              "TajD_TP_model", "TajD_TP_data", "TajD_TP_sign_match", "sign_check_required", "sign_check_pass",
              "pts", "pts_double", "corner_BP_change", "corner_TP_change", "corner_either_change",
              "numerics_run", "numerically_unreliable", "note"]


def fmt(x):
    if x is None or isinstance(x, str):
        return "NA" if x is None else x
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    return "NA" if x != x else f"{x:.6g}"


def delta_grid(n1, n2):
    i, j = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
    return i / n1 - j / n2


def corners(n1, n2):
    D = delta_grid(n1, n2)
    return {"BP": D >= CORNER - EPS, "TP": D <= -CORNER + EPS, "either": np.abs(D) >= CORNER - EPS}


def dxy_count(arr):
    """Summed between-population differences of a 2D spectrum (filled array, unfolded or folded): every
    cell's p1(1-p2) + p2(1-p1), unbiased for independent samples. Divide by callable sites for dxy."""
    n1, n2 = np.asarray(arr).shape[0] - 1, np.asarray(arr).shape[1] - 1
    p1 = (np.arange(n1 + 1) / n1)[:, None]
    p2 = (np.arange(n2 + 1) / n2)[None, :]
    return float((np.asarray(arr) * (p1 * (1 - p2) + p2 * (1 - p1))).sum())


def hudson_fst(pi1, pi2, dxy):
    return 1 - (pi1 + pi2) / 2 / dxy if dxy > 0 else float("nan")


def corner_change(e1, e2, sel):
    """Relative change of the expected count summed over a corner region (cells weighted by their count)."""
    a, b = float(e1[sel].sum()), float(e2[sel].sum())
    if a == 0:
        return 0.0 if b == 0 else float("inf")
    return abs(b - a) / abs(a)


def ratio(m, d):
    return m / d if d else float("nan")


def stats(fs, L):
    """Per-site pi and dxy, Fst and Tajima's D of a 2D spectrum (data or theta-scaled model).

    Folded spectra (checked against dadi 2.4.3 Spectrum_mod.py): marginalize() unfolds, sums and refolds;
    pi() sums n/(n-1) 2 p(1-p) and Tajima_D() uses pi() and S() (the sum of unmasked non-corner cells) over
    unmasked cells; Fst() is Weir & Cockerham's weighted ratio sum(fs*a) / sum(fs*(a+d)); dxy weights cells
    by p1(1-p2) + p2(1-p1). Every weight is unchanged when all populations are relabelled p -> 1-p, and
    fold() moves each cell's count onto its mirror cell (ambiguous cells split half and half), so all four
    are defined on folded spectra and equal their unfolded values. No statistic is NA for folded fits.
    """
    bp, tp = fs.marginalize([1]), fs.marginalize([0])
    pi_bp, pi_tp = float(bp.pi()) / L, float(tp.pi()) / L
    dxy = dxy_count(np.ma.filled(fs, 0.0)) / L
    return {"pi_BP": pi_bp, "pi_TP": pi_tp, "dxy": dxy, "Fst_WC": float(fs.Fst()),
            "Fst_Hudson": hudson_fst(pi_bp, pi_tp, dxy),
            "TajD_BP": float(bp.Tajima_D()), "TajD_TP": float(tp.Tajima_D())}


def folded_array(arr):
    return np.ma.filled(dadi.Spectrum(arr, mask_corners=True).fold(), 0.0)


def snp_free_invariant(chunk_suffix, pairs, mask_c3):
    """Fitting-mask invariant sites in 1-Mb blocks with no fitting-mask SNP, per (pair, variant), read from the
    03_spectra.py chunks as 06_gather.py reads them. None if the chunks are not there."""
    pat = re.compile(r"spectra\.\d+" + re.escape(chunk_suffix) + r"\.pkl$")
    files = sorted(f for f in glob.glob(os.path.join(OUT, "sfs/chunks/spectra.*.pkl")) if pat.search(f))
    if not files:
        return None
    snp, inv = collections.defaultdict(set), collections.defaultdict(collections.Counter)
    for path in files:
        with open(path, "rb") as fh:
            p = pickle.load(fh)
        chrom = p["chrom"]
        for var, spec_key, inv_key in (("", "spectra", "invariant"), ("_alt", "spectra_alt", "invariant_all")):
            for pid, sp in p.get(spec_key, {}).items():
                if pid in pairs:
                    snp[(pid, var)] |= {(chrom, blk) for blk, c3, region in sp
                                        if region == "none" and c3 in mask_c3}
            for (pid, blk, c3, region), n in p.get(inv_key, {}).items():
                if pid in pairs and region == "none" and c3 in mask_c3:
                    inv[(pid, var)][(chrom, blk)] += n
    return {k: sum(n for b, n in cnt.items() if b not in snp[k]) for k, cnt in inv.items()}


def spectra_rows(pair, pk, snp_free=None):
    """Section 9.2 checks on one pair's gathered spectra (sfs/<pair>.pkl). snp_free: snp_free_invariant() output."""
    rows = []
    variants = [""] + (["_alt"] if "data_alt" in pk else [])
    for v in variants:
        data = np.asarray(pk[f"data{v}"], dtype=float)
        bsum = np.asarray(pk[f"block_spec{v}"], dtype=float).sum(axis=0)
        ok = np.shape(bsum) == data.shape and np.allclose(bsum, data, rtol=RTOL, atol=ATOL)
        diff = float(np.abs(bsum - data).max()) if np.shape(bsum) == data.shape else float("nan")
        rows.append((pair, f"block_sum{v}", fmt(diff), "0 (max cell difference)", ok))
        # a block holding invariant sites but no SNP in the mask drops out of block_inv (06_gather.py keys
        # block_inv by the SNP blocks); those sites are in L_mask, so they are added back from the chunks
        extra = None if snp_free is None else snp_free.get((pair, v), 0)
        tot = data.sum() + float(np.asarray(pk[f"block_inv{v}"], dtype=float).sum()) + (extra or 0)
        L = pk[f"L_mask{v}"]
        rows.append((pair, f"total{v}", f"{tot:.1f}", f"{L:.1f} (L_mask{v})", abs(tot - L) <= RTOL * L + ATOL))
        rows.append((pair, f"invariant_in_snp_free_blocks{v}",
                      "NA (no sfs/chunks files)" if extra is None else f"{extra}",
                      "information: in L_mask, not in block_inv", "NA"))
    # Section 9.2 asks that folding the unfolded spectrum equal the directly folded spectrum. Under Deviation 1
    # the two variants hold different sites: the polarized spectrum keeps SNPs on which the outgroups agree
    # (and outgroup-callable invariant sites), the label-free spectrum every SNP, ALT-oriented. So
    # fold(polarized) cannot equal the label-free fold, and the literal check is applied within the label-free
    # variant: the folded fit target (07_fit_models.data_of) is the fold of its own unfolded spectrum. Across
    # variants only the subset relation can hold. Folding forgets orientation, and both variants project each
    # SNP from the same AN and AC, so the folded polarized spectrum must lie inside the folded label-free one
    # cell by cell, and so must its total and L.
    if "data_alt" in pk:
        direct = dadi.Spectrum(pk["data_alt"], mask_corners=True).fold()
        target = fitmod.data_of(pk, "M0f")
        same_mask = np.array_equal(np.ma.getmaskarray(direct), np.ma.getmaskarray(target))
        diff = float(np.abs(np.ma.filled(direct, 0.0) - np.ma.filled(target, 0.0)).max())
        rows.append((pair, "fold_alt_equals_fit_target", fmt(diff), "0 (max cell difference, same mask)",
                     same_mask and diff <= ATOL))
        fp, fa = folded_array(pk["data"]), np.ma.filled(direct, 0.0)
        rows.append((pair, "fold_polarized_total_le_labelfree", f"{fp.sum():.1f}", f"<= {fa.sum():.1f}",
                     fp.sum() <= fa.sum() * (1 + RTOL) + ATOL))
        rows.append((pair, "fold_polarized_cells_le_labelfree", fmt(float((fp - fa).max())), "<= 0 (max cell excess)",
                     bool(np.all(fp <= fa * (1 + RTOL) + ATOL))))
        rows.append((pair, "L_mask_polarized_le_labelfree", f"{pk['L_mask']:.1f}", f"<= {pk['L_mask_alt']:.1f}",
                     pk["L_mask"] <= pk["L_mask_alt"]))
    else:
        rows.append((pair, "fold", "not run: no label-free spectrum in this pickle", "", "NA"))
    rows.append((pair, "from_data_dict", "not run here (needs site table); see TODO",
                 "count-table spectrum == dadi.Spectrum.from_data_dict on a 1-Mb subset", "NA"))
    return rows


def model_row(pair, model, r, pk, skip_numerics):
    """Section 9.3 checks for one fitted model."""
    folded = fitmod.is_folded(model)
    v = "_alt" if folded else ""
    L = pk[f"L_mask{v}"]
    row = dict.fromkeys(MODEL_COLS)
    row.update(pair=pair, model=model, analysis="folded" if folded else "unfolded", converged=r["converged"],
               L_mask=f"{L:.0f}", numerics_run=False, numerically_unreliable=False, note="")
    bp_pop = pk["groups"][0]
    row["sign_check_required"] = bp_pop[1:] in SIGN_POPS        # wBAM, pBIK, ...: drop the sample prefix
    data = fitmod.data_of(pk, model)
    ns = data.sample_sizes
    func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(model, ns)
    popt = [r["params"][p] for p in pnames]
    row["pts"] = ",".join(str(p) for p in pts_l)
    try:
        m = func_ex(popt, ns, pts_l)
        mf = m.fold() if folded else m
        theta = dadi.Inference.optimal_sfs_scaling(mf, data)
        sm, sd = stats(theta * mf, L), stats(data, L)
        for k in ("pi_BP", "pi_TP", "dxy", "Fst_WC", "Fst_Hudson"):
            row[f"{k}_model"], row[f"{k}_data"], row[f"{k}_ratio"] = sm[k], sd[k], ratio(sm[k], sd[k])
        for side in ("BP", "TP"):
            dm, dd = sm[f"TajD_{side}"], sd[f"TajD_{side}"]
            row[f"TajD_{side}_model"], row[f"TajD_{side}_data"] = dm, dd
            row[f"TajD_{side}_sign_match"] = bool(np.isfinite(dm) and np.isfinite(dd) and np.sign(dm) == np.sign(dd))
        row["sign_check_pass"] = row["TajD_BP_sign_match"] if row["sign_check_required"] else None
    except Exception as e:                       # report, do not hide
        row["note"] = f"model statistics failed: {e}"
        row["numerically_unreliable"] = True
        return row
    if skip_numerics:
        row["note"] = "numerics skipped (--skip_numerics)"
        return row
    # section 9.3: the same parameters on twice the grid (pts = [2n+20, 2n+40, 2n+60] with model_setup's
    # [n+10, n+20, n+30]). E is rebuilt exactly as 09_tail.py would build it on the finer grid, theta included,
    # in the unfolded shape that 09_tail.py sums its tail over.
    pts2 = [2 * p for p in pts_l]
    row["pts_double"] = ",".join(str(p) for p in pts2)
    try:
        t0 = time.time()
        m2 = func_ex(popt, ns, pts2)
        th2 = dadi.Inference.optimal_sfs_scaling(m2.fold() if folded else m2, data)
        e1, e2 = np.ma.filled(theta * m, 0.0), np.ma.filled(th2 * m2, 0.0)
        sel = corners(*ns)
        keys = ("either",) if folded else ("BP", "TP")
        for k in keys:
            row[f"corner_{k}_change"] = corner_change(e1, e2, sel[k])
        row["numerics_run"] = True
        row["numerically_unreliable"] = any(row[f"corner_{k}_change"] > MAX_CHANGE for k in keys)
        row["note"] = f"doubled grid {time.time() - t0:.0f} s"
    except Exception as e:                       # an unverifiable grid counts against the model
        row["note"] = f"doubled grid failed: {e}"
        row["numerically_unreliable"] = True
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="")
    ap.add_argument("--chunk_suffix", default=None, help="suffix of the 03_spectra.py chunks; defaults to --suffix")
    ap.add_argument("--pairs", default="A1w,A2w,A3w,A4w,C1w,C2w")
    ap.add_argument("--skip_numerics", action="store_true",
                    help="skip the doubled-grid corner check (slow); numerically_unreliable is then False")
    a = ap.parse_args()

    pairs = a.pairs.split(",")
    bf = os.path.join(OUT, f"fits/best_fits{a.suffix}.tsv")
    for path in [bf] + [os.path.join(OUT, f"sfs/{p}{a.suffix}.pkl") for p in pairs]:
        if not os.path.exists(path):
            raise SystemExit(f"missing input {path}")
    best = collections.defaultdict(dict)
    with open(bf) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["params"] = {k: float(v) for k, v in (x.split("=") for x in r["params"].split(";"))}
            best[r["pair"]][r["model"]] = r

    with open(os.path.join(OUT, f"sfs/{pairs[0]}{a.suffix}.pkl"), "rb") as fh:
        mask_c3 = set(pickle.load(fh)["mask_c3"])
    snp_free = snp_free_invariant(a.suffix if a.chunk_suffix is None else a.chunk_suffix, set(pairs), mask_c3)
    srows, mrows = [], []
    for pair in pairs:
        with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.pkl"), "rb") as fh:
            pk = pickle.load(fh)
        if set(pk["mask_c3"]) != mask_c3:
            raise SystemExit(f"{pair}: mask_c3 {pk['mask_c3']} differs from {sorted(mask_c3)}")
        sr = spectra_rows(pair, pk, snp_free)
        srows += sr
        checked = [x for x in sr if x[4] != "NA"]
        failed = [x[1] for x in checked if not x[4]]
        mr = []
        for model, r in sorted(best.get(pair, {}).items()):
            if fitmod.is_folded(model) and "data_alt" not in pk:     # runs gathered before the label-free spectra
                continue
            mr.append(model_row(pair, model, r, pk, a.skip_numerics))
        mrows += mr
        req = [x for x in mr if x["sign_check_required"] and x["sign_check_pass"] is not None]
        bad = [x["model"] for x in mr if x["numerically_unreliable"]]
        print(f"{pair}: spectra {len(checked) - len(failed)}/{len(checked)} pass"
              + (f" (FAIL {','.join(failed)})" if failed else "")
              + f"; {len(mr)} models" + ("" if pair in best else " (no fits)")
              + (f"; TajD BP sign {sum(x['sign_check_pass'] for x in req)}/{len(req)} match" if req else "")
              + ("; numerics skipped" if a.skip_numerics else
                 f"; numerically unreliable: {','.join(bad) if bad else 'none'}"), flush=True)

    os.makedirs(os.path.join(OUT, "checks"), exist_ok=True)
    with open(os.path.join(OUT, f"checks/spectra_checks{a.suffix}.tsv"), "w") as o:
        o.write("pair\tcheck\tvalue\texpected\tpass\n")
        for r in srows:
            o.write("\t".join(fmt(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"checks/model_checks{a.suffix}.tsv"), "w") as o:
        o.write("\t".join(MODEL_COLS) + "\n")
        for r in mrows:
            o.write("\t".join(fmt(r[c]) for c in MODEL_COLS) + "\n")
    print(f"wrote checks/spectra_checks{a.suffix}.tsv ({len(srows)} rows), "
          f"checks/model_checks{a.suffix}.tsv ({len(mrows)} rows)")


if __name__ == "__main__":
    main()
