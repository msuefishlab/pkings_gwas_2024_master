#!/usr/bin/env python3
"""09_tail.py
How many SNPs as BP-biased as each index SNP does drift make, under the fitted histories?
DESIGN.md sections 4 and 6 and Deviation 1, applied exactly as pre-specified.

Two analyses per pair:
  unfolded  models M0-M3 fitted to the outgroup-polarized spectrum (with a misidentification parameter);
            tail R = projected cells with i/n_BP - j/n_TP >= delta_idx (derived-allele delta)
  folded    models M0f-M3f fitted to the label-free spectrum of every SNP; tail R = cells with
            |i/n_BP - j/n_TP| >= |delta_alt| in either direction, so its E counts both directions and
            can only exceed the directional E (drift-generous)
For each model:
  exp          theta-scaled model spectrum over the fitting mask (unfolded shape; folded fits use the
               demographic model without misidentification)
  ladder       observed and expected counts at delta >= {0.5 .. 1.0}
  admissible   converged, p_misid <= 0.05 (unfolded only), and at the most extreme ladder threshold with
               E >= 10, 1/3 <= O/E <= 3
  tau          each branch's total drift since the split
Headline (Deviation 1): the unfolded analysis if any unfolded model is admissible for the pair; otherwise
the folded analysis. Both are reported.
For each unit (pair x index SNP) and analysis:
  E_genome     E over R scaled from the fitting mask to every callable site (L_scan / L_mask); for a converged
               model, the largest E among its best three starts (Deviation 2: the starts are statistically
               indistinguishable, so the most drift-generous is used)
  E_hetNe      the same with local Ne scaled by f_w per 1-Mb block (T/f, m*f, theta*f)
  E_robust     max(upper 97.5% bootstrap E under the best model, point E under every admissible model,
               E_hetNe under the best model)
  tier         <= 0.005, <= 0.05, <= 1, <= 20, > 20, with the section 6.3 overrides
  K_needed     L_scan / E_robust
  f_star       largest grid f at which the peak region (~500 kb) alone expects >= 1 such site; 0 = none

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
EPS = 1e-9


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


def rescale(pnames, popt, f):
    """Local Ne x f: times in 2N units divide by f, migration (2Nm) multiplies by f."""
    p = dict(zip(pnames, popt))
    p["T"] = p["T"] / f
    for k in ("m", "m_BP_from_TP", "m_TP_from_BP"):
        if k in p:
            p[k] = p[k] * f
    return [p[k] for k in pnames]


def taus(model, p):
    T = p["T"]
    if model.rstrip("f") == "M3":
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
            r["top3"] = ([{k: float(v) for k, v in (x.split("=") for x in t.split(";"))}
                          for t in r["top3_params"].split("|")] if r.get("top3_params") else [r["params"]])
            best[r["pair"]][r["model"]] = r
    boot = collections.defaultdict(list)
    for path in glob.glob(os.path.join(OUT, f"fits/*{a.suffix}.M*.boot*.seed*.tsv")):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                boot[(r["pair"], r["model"])].append(r)
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
        n1, n2 = pk["n"]
        i_idx, j_idx = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
        D = i_idx / n1 - j_idx / n2

        def variant(folded):
            v = "_alt" if folded else ""
            lw = pk[f"block_spec{v}"].sum(axis=(1, 2)) + pk[f"block_inv{v}"]
            fw = np.array([pk["f_blocks"].get(k, 1.0) for k in pk[f"block_keys{v}"]])
            return dict(obs=np.asarray(pk[f"data{v}"]), L_mask=pk[f"L_mask{v}"], L_scan=pk[f"L_scan{v}"],
                        region_spec=pk[f"region_spec{v}"], lw=lw, fw=fw)
        VAR = {False: variant(False), True: variant(True) if "data_alt" in pk else None}

        def sel(folded, thr):
            return (np.abs(D) >= thr - EPS) if folded else (D >= thr - EPS)

        info = {}
        for model, r in sorted(best[pair].items()):
            folded = fitmod.is_folded(model)
            if VAR[folded] is None:      # runs gathered before the label-free spectra existed
                continue
            V = VAR[folded]
            data = fitmod.data_of(pk, model)
            func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(model, data.sample_sizes)
            popt = [r["params"][p] for p in pnames]
            m = func_ex(popt, data.sample_sizes, pts_l)
            theta = dadi.Inference.optimal_sfs_scaling(m.fold() if folded else m, data)
            exp = np.ma.filled(theta * m, 0.0)                    # unfolded shape, expected counts
            exp_ridge = [exp]                                      # Deviation 2: best three starts
            if r["converged"]:
                for t in r["top3"][1:]:
                    mt = func_ex([t[p_] for p_ in pnames], data.sample_sizes, pts_l)
                    th = dadi.Inference.optimal_sfs_scaling(mt.fold() if folded else mt, data)
                    exp_ridge.append(np.ma.filled(th * mt, 0.0))
            ladder = []
            for thr in LADDER:
                checks = [("either" if folded else "BP_derived", sel(folded, thr))]
                if not folded:
                    checks.append(("TP_derived", D <= -thr + EPS))
                for direction, s_ in checks:
                    O, E = V["obs"][s_].sum(), exp[s_].sum()
                    ladder.append((thr, direction, O, E))
                    lrows.append((pair, model, thr, direction, O, E, O / E if E else float("nan")))
            gate = [(thr, O, E) for thr, d, O, E in ladder if d in ("BP_derived", "either") and E >= 10]
            thr_g, O_g, E_g = gate[-1] if gate else (float("nan"), float("nan"), float("nan"))
            adequate = bool(gate) and (1 / 3) <= O_g / E_g <= 3
            p = dict(zip(pnames, popt))
            misid = p.get("p_misid", float("nan"))
            admissible = r["converged"] and adequate and (folded or misid <= 0.05)
            tb, tt = taus(model, p)
            q_grid = {}
            for f in F_GRID:
                mf = func_ex(rescale(pnames, popt, f), data.sample_sizes, pts_l)
                q_grid[f] = np.ma.filled(theta * f * mf, 0.0) / V["L_mask"]
            info[model] = dict(r=r, folded=folded, popt=popt, pnames=pnames, func_ex=func_ex, pts_l=pts_l,
                               ns=data.sample_sizes, exp=exp, exp_ridge=exp_ridge, theta=theta,
                               admissible=admissible, q_grid=q_grid)
            n_anc = [theta / (4 * mu * V["L_mask"]) for mu in MU]
            bp95 = D >= 0.95 - EPS
            tp95 = D <= -0.95 + EPS
            mrows.append((pair, model, "folded" if folded else "unfolded", r["LL"], r["converged"],
                          f"{misid:.4f}", thr_g, f"{O_g:.1f}", f"{E_g:.1f}", adequate, admissible, f"{theta:.1f}",
                          f"{p['T']:.5f}", f"{tb:.5f}", f"{tt:.5f}", f"{1 - math.exp(-tb):.4f}",
                          f"{1 - math.exp(-tt):.4f}", f"{n_anc[0]:.0f}", f"{n_anc[1]:.0f}",
                          f"{V['obs'][bp95].sum():.1f}", f"{exp[bp95].sum():.1f}",
                          f"{V['obs'][tp95].sum():.1f}", f"{exp[tp95].sum():.1f}",
                          ";".join(f"{k}={v:.5g}" for k, v in p.items())))

        adm = {False: [m for m, v in info.items() if v["admissible"] and not v["folded"]],
               True: [m for m, v in info.items() if v["admissible"] and v["folded"]]}
        headline = False if adm[False] else True

        def analyse(u, folded):
            V = VAR[folded]
            key = "delta_alt" if folded else "delta"
            if u[key] in ("NA", "nan"):
                return None
            thr = abs(float(u[key])) if folded else float(u[key])
            R = sel(folded, thr)
            own_region = f"peak{u['peak']}"
            O_mask = V["obs"][R].sum()
            O_genome = sum(np.asarray(arr)[R].sum() for k, arr in V["region_spec"].items() if k != own_region)
            ratio = V["L_scan"] / V["L_mask"]
            models = [m for m in info if info[m]["folded"] == folded]
            E_by = {m: max(e[R].sum() for e in info[m]["exp_ridge"]) * ratio for m in models}
            res = dict(thr=thr, O_mask=O_mask, O_genome=O_genome, E_by=E_by, best=None, E_point=float("nan"),
                       E_het=float("nan"), E_boot=float("nan"), f_star=float("nan"), E_robust=float("nan"),
                       tier="model-dependent (no admissible model): report O")
            if not adm[folded]:
                return res
            bm = max(adm[folded], key=lambda m: float(info[m]["r"]["LL"]))
            v = info[bm]
            res["best"], res["E_point"] = bm, E_by[bm]
            qf = np.array([v["q_grid"][f][R].sum() for f in F_GRID])
            q_w = np.interp(np.clip(V["fw"], F_GRID[0], F_GRID[-1]), F_GRID, qf)
            res["E_het"] = float((V["lw"] * q_w).sum() * ratio)
            need = qf * V["L_scan"] * 500_000 / GENOME
            ok_f = [f for f, nd in zip(F_GRID, need) if nd >= 1]
            res["f_star"] = float(max(ok_f)) if ok_f else 0.0
            Eb = []
            for br in boot.get((pair, bm), []):
                pb = [float(br[p]) for p in v["pnames"]]
                mb = v["func_ex"](pb, v["ns"], v["pts_l"])
                Eb.append(float(br["theta"]) * np.ma.filled(mb, 0.0)[R].sum() * ratio)
            if len(Eb) >= 20:
                res["E_boot"] = float(np.quantile(Eb, 0.975))
            else:
                print(f"WARNING {pair} peak {u['peak']} {'folded' if folded else 'unfolded'}: {len(Eb)} bootstrap "
                      f"refits for best model {bm}; E_robust lacks the bootstrap bound", flush=True)
            cands = [E_by[m] for m in adm[folded]] + [x for x in (res["E_boot"], res["E_het"]) if x == x]
            res["E_robust"] = max(cands)
            tier = tier_of(res["E_robust"])
            if res["E_robust"] <= 1 and O_genome > max(5, 3 * res["E_point"]):
                tier = f"model fails at the extremes: read from O ({O_genome:.0f})"
            if len({tier_of(E_by[m]) for m in adm[folded]}) > 1:
                tier += " [model-dependent]"
            res["tier"] = tier
            return res

        for u in units.get(pair, []):
            f_peak = pk["f_blocks"].get((u["chrom"], int(u["pos"]) // 1_000_000), float("nan"))
            for folded in (False, True):
                if VAR[folded] is None:
                    continue
                res = analyse(u, folded)
                if res is None:
                    continue
                L_scan = VAR[folded]["L_scan"]
                Er = res["E_robust"]
                urows.append((pair, u["peak"], u["own"], "folded" if folded else "unfolded",
                              folded == headline, u["f_bp"], u["f_tp"], f"{res['thr']:.4f}",
                              f"{res['O_mask']:.2f}", f"{res['O_genome']:.2f}", res["best"] or "none",
                              f"{res['E_point']:.4g}", ";".join(f"{m}={e:.4g}" for m, e in res["E_by"].items()),
                              f"{res['E_het']:.4g}", f"{res['E_boot']:.4g}", f"{Er:.4g}", res["tier"],
                              f"{L_scan / Er:.4g}" if Er == Er and Er > 0 else "inf",
                              f"{f_peak:.3f}", f"{res['f_star']:.3f}"))

    os.makedirs(os.path.join(OUT, "tail"), exist_ok=True)
    with open(os.path.join(OUT, f"tail/models{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tanalysis\tLL\tconverged\tp_misid\tgate_threshold\tgate_O\tgate_E\tadequate\t"
                "admissible\ttheta\tT_2Nanc\ttau_BP\ttau_TP\tF_BP\tF_TP\tN_anc_mu4e-9\tN_anc_mu1e-8\t"
                "obs_BPder95\texp_BPder95\tobs_TPder95\texp_TPder95\tparams\n")
        for r in mrows:
            o.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/ladder{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tthreshold\tdirection\tO\tE\tO_over_E\n")
        for r in lrows:
            o.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/unit_results{a.suffix}.tsv"), "w") as o:
        o.write("pair\tpeak\town\tanalysis\theadline\tf_bp\tf_tp\ttail_threshold\tO_mask\tO_genome\tbest_model\t"
                "E_genome_best\tE_by_model\tE_hetNe\tE_boot97.5\tE_robust\ttier\tK_needed\tf_peak\tf_star\n")
        for r in urows:
            o.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote tail/models{a.suffix}.tsv, ladder, unit_results ({len(urows)} rows)")


if __name__ == "__main__":
    main()
