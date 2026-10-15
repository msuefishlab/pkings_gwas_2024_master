#!/usr/bin/env python3
"""09_tail.py
How many SNPs as BP-biased as each index SNP does drift make, under the fitted histories?
DESIGN.md sections 4 and 6, applied exactly as pre-specified.

For each pair and each fitted model:
  exp          theta-scaled, misidentification-wrapped model spectrum over the fitting mask
  ladder       observed and expected counts at delta >= {0.5 .. 1.0} (and the mirror for TP-derived)
  admissible   converged, p_misid <= 0.05, and at the most extreme ladder threshold with E >= 10,
               0.33 <= O/E <= 3
  tau          each branch's total drift since the split, F = 1 - exp(-tau)
For each unit (pair x index SNP):
  R            projected cells with i/n_BP - j/n_TP >= delta_idx (raw delta at the index SNP)
  E_mask       sum of exp over R; E_genome = E_mask * L_scan / L_mask
  E_hetNe      the same with local Ne scaled by f_w per 1-Mb block (T/f, m*f, theta*f)
  E_robust     max(upper 97.5% bootstrap E_genome under the best model, point E_genome under every
               admissible model, E_hetNe under the best model)
  tier         <= 0.005, <= 0.05, <= 1, <= 20, > 20, with the section 6.3 overrides
  f_star       largest grid f at which the peak region (~500 kb) alone expects >= 1 such site; 0 = none
  K_needed     L_scan / E_robust

Writes tail/models<suffix>.tsv, tail/ladder<suffix>.tsv, tail/unit_results<suffix>.tsv.
"""
import argparse
import collections
import csv
import glob
import importlib.util
import math
import os
import pickle

import numpy as np
import dadi

OUT = os.path.join(os.getcwd(), os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("fitmod", os.path.join(HERE, "07_fit_models.py"))
fitmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitmod)

LADDER = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
F_GRID = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.15, 1.3, 1.5, 2.0])
MU = (4e-9, 1e-8)
GENOME = 888_337_685


def tier_of(e):
    if e <= 0.005:
        return "drift implausible (all 11 units)"
    if e <= 0.05:
        return "drift implausible (this unit)"
    if e <= 1:
        return "unusual under drift"
    if e <= 20:
        return "drift makes a few"
    return "drift readily makes this"


def rescale(model, pnames, popt, f):
    """Local Ne x f: times in 2N units divide by f, migration (2Nm) multiplies by f."""
    p = dict(zip(pnames, popt))
    p["T"] = p["T"] / f
    for k in ("m", "m_BP_from_TP", "m_TP_from_BP"):
        if k in p:
            p[k] = p[k] * f
    return [p[k] for k in pnames]


def taus(model, p):
    T = p["T"]
    if model == "M3":
        def grow(start, end):
            return T / start if abs(end - start) < 1e-12 else T * (1 / start - 1 / end) / math.log(end / start)
        return grow(p["s_BP"], p["nu_BP"]), grow(1 - p["s_BP"], p["nu_TP"])
    return T / p["nu_BP"], T / p["nu_TP"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="")
    ap.add_argument("--pairs", default="A1w,A2w,A3w,A4w,C1w,C2w")
    a = ap.parse_args()

    best = collections.defaultdict(dict)
    with open(os.path.join(OUT, f"fits/best_fits{a.suffix}.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["params"] = {k: float(v) for k, v in (x.split("=") for x in r["params"].split(";"))}
            r["converged"] = r["converged"] == "True"
            best[r["pair"]][r["model"]] = r
    boot = collections.defaultdict(list)
    for path in glob.glob(os.path.join(OUT, f"fits/*{a.suffix}.M?.boot*.seed*.tsv")):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                boot[(r["pair"], r["model"], r["target"])].append(r)
    units = collections.defaultdict(list)
    with open(os.path.join(OUT, f"tail/index_units{a.suffix}.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            units[r["pair_id"]].append(r)

    mrows, lrows, urows = [], [], []
    for pair in a.pairs.split(","):
        if pair not in best:
            print("no fits for", pair)
            continue
        with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.pkl"), "rb") as fh:
            pk = pickle.load(fh)
        data = dadi.Spectrum(pk["data"], mask_corners=True)
        n1, n2 = pk["n"]
        ns = data.sample_sizes
        i_idx, j_idx = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
        D = i_idx / n1 - j_idx / n2
        L_mask, L_scan = pk["L_mask"], pk["L_scan"]
        lw = pk["block_spec"].sum(axis=(1, 2)) + pk["block_inv"]
        fw = np.array([pk["f_blocks"].get(k, 1.0) for k in pk["block_keys"]])

        model_info = {}
        for model, r in sorted(best[pair].items()):
            func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(model, ns)
            popt = [r["params"][p] for p in pnames]
            m = func_ex(popt, ns, pts_l)
            theta = dadi.Inference.optimal_sfs_scaling(m, data)
            exp = np.ma.filled(theta * m, 0.0)
            obs = np.ma.filled(data, 0.0)
            ladder = []
            for thr in LADDER:
                for direction, sel in (("BP_derived", D >= thr - 1e-9), ("TP_derived", D <= -thr + 1e-9)):
                    O, E = obs[sel].sum(), exp[sel].sum()
                    ladder.append((thr, direction, O, E))
                    lrows.append((pair, model, thr, direction, O, E, O / E if E else float("nan")))
            gate = [(thr, O, E) for thr, d, O, E in ladder if d == "BP_derived" and E >= 10]
            thr_g, O_g, E_g = gate[-1] if gate else (float("nan"), float("nan"), float("nan"))
            adequate = bool(gate) and (1 / 3) <= O_g / E_g <= 3
            p = dict(zip(pnames, popt))
            tb, tt = taus(model, p)
            admissible = r["converged"] and p["p_misid"] <= 0.05 and adequate
            n_anc = [theta / (4 * mu * L_mask) for mu in MU]
            # q(f) over the f grid, for E_hetNe and f*
            q_grid = {}
            for f in F_GRID:
                mf = func_ex(rescale(model, pnames, popt, f), ns, pts_l)
                q_grid[f] = np.ma.filled(theta * f * mf, 0.0) / L_mask
            model_info[model] = dict(r=r, popt=popt, pnames=pnames, func_ex=func_ex, pts_l=pts_l, exp=exp,
                                     theta=theta, admissible=admissible, tau=(tb, tt), q_grid=q_grid)
            obs_bp = obs[D >= 0.95 - 1e-9].sum(); obs_tp = obs[D <= -0.95 + 1e-9].sum()
            exp_bp = exp[D >= 0.95 - 1e-9].sum(); exp_tp = exp[D <= -0.95 + 1e-9].sum()
            mrows.append((pair, model, f"{r['LL']}", r["converged"], f"{p['p_misid']:.4f}", f"{thr_g}",
                          f"{O_g:.1f}", f"{E_g:.1f}", adequate, admissible, f"{theta:.1f}",
                          f"{p['T']:.5f}", f"{tb:.5f}", f"{tt:.5f}", f"{1 - math.exp(-tb):.4f}", f"{1 - math.exp(-tt):.4f}",
                          f"{n_anc[0]:.0f}", f"{n_anc[1]:.0f}", f"{obs_bp:.1f}", f"{exp_bp:.1f}", f"{obs_tp:.1f}", f"{exp_tp:.1f}",
                          ";".join(f"{k}={v:.5g}" for k, v in p.items())))

        adm = {m: v for m, v in model_info.items() if v["admissible"]}
        best_model = max(adm, key=lambda m: float(adm[m]["r"]["LL"])) if adm else None

        for u in units.get(pair, []):
            if u["delta"] == "NA":
                continue
            delta = float(u["delta"])
            R = D >= delta - 1e-9
            own_region = f"peak{u['peak']}"
            O_mask = np.ma.filled(data, 0.0)[R].sum()
            O_genome = sum(np.asarray(v)[R].sum() for k, v in pk["region_spec"].items() if k != own_region)
            E_by_model = {}
            for model, v in model_info.items():
                E_by_model[model] = v["exp"][R].sum() * L_scan / L_mask
            E_het = float("nan"); E_boot_hi = float("nan"); f_star = float("nan"); E_point = float("nan")
            if best_model:
                v = model_info[best_model]
                E_point = E_by_model[best_model]
                qf = np.array([v["q_grid"][f][R].sum() for f in F_GRID])
                q_w = np.interp(np.clip(fw, F_GRID[0], F_GRID[-1]), F_GRID, qf)
                E_het = float((lw * q_w).sum() * L_scan / L_mask)
                # f*: local Ne at which the peak region alone expects one such site
                L_peak = L_scan * 500_000 / GENOME
                need = qf * L_peak
                # the largest grid f at which the peak region alone expects >= 1 such site (0 if none)
                ok_f = [f for f, nd in zip(F_GRID, need) if nd >= 1]
                f_star = float(max(ok_f)) if ok_f else 0.0
                Eb = []
                for br in [x for k, xs in boot.items() if k[0] == pair and k[1] == best_model for x in xs]:
                    pb = [float(br[p]) for p in v["pnames"]]
                    mb = v["func_ex"](pb, ns, v["pts_l"])
                    Eb.append(float(br["theta"]) * np.ma.filled(mb, 0.0)[R].sum() * L_scan / L_mask)
                if len(Eb) >= 20:
                    E_boot_hi = float(np.quantile(Eb, 0.975))
            cands = [E_by_model[m] for m in adm] + [x for x in (E_boot_hi, E_het) if x == x]
            E_robust = max(cands) if cands else float("nan")
            if not adm:
                tier = "model-dependent (no admissible model): report O"
            else:
                tier = tier_of(E_robust)
                if E_robust <= 1 and O_genome > max(5, 3 * E_point):
                    tier = f"model fails at the extremes: read from O ({O_genome:.0f})"
                tiers_all = {tier_of(E_by_model[m]) for m in adm}
                if len(tiers_all) > 1:
                    tier += " [model-dependent]"
            f_peak = pk["f_blocks"].get((u["chrom"], int(u["pos"]) // 1_000_000), float("nan"))
            urows.append((pair, u["peak"], u["own"], u["f_bp"], u["f_tp"], f"{delta:.4f}", f"{O_mask:.2f}",
                          f"{O_genome:.2f}", best_model or "none", f"{E_point:.4g}",
                          ";".join(f"{m}={e:.4g}" for m, e in E_by_model.items()), f"{E_het:.4g}", f"{E_boot_hi:.4g}",
                          f"{E_robust:.4g}", tier, f"{L_scan / E_robust:.4g}" if E_robust == E_robust and E_robust > 0 else "inf",
                          f"{f_peak:.3f}", f"{f_star:.3f}"))

    os.makedirs(os.path.join(OUT, "tail"), exist_ok=True)
    with open(os.path.join(OUT, f"tail/models{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tLL\tconverged\tp_misid\tgate_threshold\tgate_O\tgate_E\tadequate\tadmissible\ttheta\t"
                "T_2Nanc\ttau_BP\ttau_TP\tF_BP\tF_TP\tN_anc_mu4e-9\tN_anc_mu1e-8\tobs_BPder95\texp_BPder95\t"
                "obs_TPder95\texp_TPder95\tparams\n")
        for r in mrows:
            o.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/ladder{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tthreshold\tdirection\tO\tE\tO_over_E\n")
        for r in lrows:
            o.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/unit_results{a.suffix}.tsv"), "w") as o:
        o.write("pair\tpeak\town\tf_bp\tf_tp\tdelta\tO_mask\tO_genome\tbest_model\tE_genome_best\tE_by_model\t"
                "E_hetNe\tE_boot97.5\tE_robust\ttier\tK_needed\tf_peak\tf_star\n")
        for r in urows:
            o.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote tail/models{a.suffix}.tsv, ladder, unit_results ({len(urows)} units)")


if __name__ == "__main__":
    main()
