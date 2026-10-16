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
Headline per unit (Addenda 1a, 1b): the unfolded analysis if the pair has an admissible unfolded model and the
index SNP is polarized; otherwise the folded analysis. Both are reported.
For each unit (pair x index SNP) and analysis:
  E_genome     E over R scaled from the fitting mask to every callable site (L_scan / L_mask); for a converged
               model, the largest E among its best three starts (Deviation 2: the starts are statistically
               indistinguishable, so the most drift-generous is used)
  E_hetNe      the same with local Ne scaled by f_w per 1-Mb block (T/f, m*f, theta*f), f clipped to 0.2-1.5,
               each block weighted by its callable sites over the whole scan (06b_block_scan.py); the same on
               the wider f grid 0.05-2.0 is reported as information
  E_grid2      E under the best model with the grid doubled (pts -> 2 pts; Addendum 2b)
  E_robust     max(upper 97.5% bootstrap E under the best model, point E under every admissible model,
               E_hetNe and E_grid2 under the best model)
  tier         <= 0.005, <= 0.05, <= 1, <= 20, > 20, with the section 6.3 overrides
  K_needed     L_scan / E_genome under the best model (section 4); K_needed_robust uses E_robust
  f_star       largest grid f at which the peak region (~500 kb) alone expects >= 1 such site; 0 = none

Writes tail/models<suffix>.tsv (with bootstrap intervals for T/2N and tau), tail/ladder<suffix>.tsv,
tail/unit_results<suffix>.tsv and tail/drift_exposure<suffix>.tsv (section 6.4 inputs).
"""
import argparse
import collections
import csv
import glob
import importlib.util
import math
import multiprocessing
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
F_DESIGN = (0.2, 1.5)    # DESIGN.md section 4: E_hetNe on f in 0.2-1.5; the wider grid is reported as information
MU = (4e-9, 1e-8)
GENOME = 888_337_685
# index_units.tsv stores delta to 4 decimals, so a delta on the projected grid can be rounded up past its own
# cell. Projected grid values are multiples of 1/lcm(n_BP, n_TP) >= 1/1380 apart, so 1e-4 keeps that cell in R
# without admitting the next one.
EPS = 1e-4


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


def _do_pair(pair):
    """One pair's rows of the six output tables. Runs in a worker process and reads the inputs that main()
    sets as module globals (a, best, boot, sims, unreliable, units)."""
    mrows, lrows, urows, drows, rrows, crows = [], [], [], [], [], []
    if pair not in best:
        print("no fits for", pair)
        return mrows, lrows, urows, drows, rrows, crows
    with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.pkl"), "rb") as fh:
        pk = pickle.load(fh)
    n1, n2 = pk["n"]
    i_idx, j_idx = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
    D = i_idx / n1 - j_idx / n2

    # E_hetNe weights: callable sites per 1-Mb block over the whole scan (06b_block_scan.py), not the fitting
    # mask, which leaves out the exon-dense, low-f blocks where background selection is strongest
    scan = collections.defaultdict(dict)
    with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.block_scan.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            scan[""][(r["chrom"], int(r["block"]))] = float(r["L_scan"])
            scan["_alt"][(r["chrom"], int(r["block"]))] = float(r["L_scan_alt"])

    def variant(folded):
        v = "_alt" if folded else ""
        keys = sorted(scan[v])
        lw = np.array([scan[v][k] for k in keys])
        if abs(lw.sum() - pk[f"L_scan{v}"]) > 1e-6 * pk[f"L_scan{v}"]:
            raise SystemExit(f"{pair}{v}: block_scan sums to {lw.sum():.0f}, L_scan is {pk[f'L_scan{v}']:.0f}")
        fw = np.array([pk["f_blocks"].get(k, 1.0) for k in keys])
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
        info[model] = dict(r=r, folded=folded, popt=popt, pnames=pnames, func_ex=func_ex, pts_l=pts_l,
                           data=data, ns=data.sample_sizes, exp=exp, exp_ridge=exp_ridge, theta=theta,
                           admissible=admissible)
        n_anc = [theta / (4 * mu * V["L_mask"]) for mu in MU]
        # asymmetry (section 4): derived-allele counts, so unfolded models only; a folded model's spectrum is
        # ALT-oriented and has no BP-derived or TP-derived side
        bp95 = D >= 0.95 - EPS
        tp95 = D <= -0.95 + EPS
        asym = (["NA"] * 4 if folded else
                [f"{V['obs'][bp95].sum():.1f}", f"{exp[bp95].sum():.1f}",
                 f"{V['obs'][tp95].sum():.1f}", f"{exp[tp95].sum():.1f}"])
        # section 5: bootstrap refits (best model and M3) give the intervals for T/2N and tau
        ci = {"T": "NA", "tau_BP": "NA", "tau_TP": "NA"}
        brs = boot.get((pair, model), [])
        if len(brs) >= 20:
            bt = np.array([[float(b["T"]), *taus(model, {k: float(b[k]) for k in pnames})] for b in brs])
            for col, name in enumerate(("T", "tau_BP", "tau_TP")):
                lo, hi = np.quantile(bt[:, col], [0.025, 0.975])
                ci[name] = f"{lo:.5f}-{hi:.5f}"
        mrows.append((pair, model, "folded" if folded else "unfolded", r["LL"], r["converged"],
                      f"{misid:.4f}", thr_g, f"{O_g:.1f}", f"{E_g:.1f}", adequate, admissible, f"{theta:.1f}",
                      f"{p['T']:.5f}", f"{tb:.5f}", f"{tt:.5f}", f"{1 - math.exp(-tb):.4f}",
                      f"{1 - math.exp(-tt):.4f}", f"{n_anc[0]:.0f}", f"{n_anc[1]:.0f}", *asym,
                      len(brs), ci["T"], ci["tau_BP"], ci["tau_TP"],
                      ";".join(f"{k}={v:.5g}" for k, v in p.items())))
        info[model].update(tau_BP=tb, tau_TP=tt, asym=asym)
        # section 9.3 parameter recovery: is each generating value inside the 95% range of the Poisson refits?
        srs = sims.get((pair, model), [])
        if len(srs) >= 20:
            for k in pnames:
                est = np.array([float(s_[k]) for s_ in srs])
                lo, med, hi = np.quantile(est, [0.025, 0.5, 0.975])
                rrows.append((pair, model, k, "", len(srs), f"{p[k]:.6g}", f"{lo:.6g}", f"{med:.6g}",
                              f"{hi:.6g}", lo <= p[k] <= hi))

    adm = {False: [m for m, v in info.items() if v["admissible"] and not v["folded"]],
           True: [m for m, v in info.items() if v["admissible"] and v["folded"]]}

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
                   E_het=float("nan"), E_het_wide=float("nan"), E_boot=float("nan"), E_grid2=float("nan"),
                   f_star=float("nan"),
                   E_robust=float("nan"), tier="model-dependent (no admissible model): report O")
        if not adm[folded]:
            return res
        bm = max(adm[folded], key=lambda m: float(info[m]["r"]["LL"]))
        v = info[bm]
        res["best"], res["E_point"] = bm, E_by[bm]
        if "q_grid" not in v:        # per-site spectrum at local Ne x f; built once, for best models only
            v["q_grid"] = {f: np.ma.filled(v["theta"] * f * v["func_ex"](rescale(v["pnames"], v["popt"], f),
                                                                         v["ns"], v["pts_l"]), 0.0) / V["L_mask"]
                           for f in F_GRID}
        qf = np.array([v["q_grid"][f][R].sum() for f in F_GRID])
        # sum over every scanned block of its callable sites x the per-site tail probability at its local Ne
        res["E_het"] = float((V["lw"] * np.interp(np.clip(V["fw"], *F_DESIGN), F_GRID, qf)).sum())
        res["E_het_wide"] = float((V["lw"] * np.interp(np.clip(V["fw"], F_GRID[0], F_GRID[-1]), F_GRID, qf)).sum())
        need = qf * V["L_scan"] * 500_000 / GENOME
        ok_f = [f for f, nd in zip(F_GRID, need) if nd >= 1]
        res["f_star"] = float(max(ok_f)) if ok_f else 0.0
        # Addendum 2b: corner cells can move > 5% when the grid doubles (section 9.3), so E on twice the grid
        # (theta refitted there) is a further candidate for E_robust; computed once per best model
        if "exp_grid2" not in v:
            m2 = v["func_ex"](v["popt"], v["ns"], [2 * p_ for p_ in v["pts_l"]])
            th2 = dadi.Inference.optimal_sfs_scaling(m2.fold() if folded else m2, v["data"])
            v["exp_grid2"] = np.ma.filled(th2 * m2, 0.0)
        res["E_grid2"] = float(v["exp_grid2"][R].sum() * ratio)
        # theta-scaled spectra of the bootstrap and Poisson refits, built once per best model and reused for
        # every unit (each evaluation costs seconds at 60 x 60)
        def refit_spectra(rows):
            return [float(s_["theta"]) * np.ma.filled(v["func_ex"]([float(s_[p_]) for p_ in v["pnames"]],
                                                                    v["ns"], v["pts_l"]), 0.0) for s_ in rows]
        if "boot_exp" not in v:
            v["boot_exp"] = refit_spectra(boot.get((pair, bm), []))
        Eb = [e[R].sum() * ratio for e in v["boot_exp"]]
        if len(Eb) >= 20:
            res["E_boot"] = float(np.quantile(Eb, 0.975))
        else:
            print(f"WARNING {pair} peak {u['peak']} {'folded' if folded else 'unfolded'}: {len(Eb)} bootstrap "
                  f"refits for best model {bm}; E_robust lacks the bootstrap bound", flush=True)
        srs = sims.get((pair, bm), [])
        if len(srs) >= 20:              # section 9.3: E recovered from the Poisson refits
            e_true = v["exp"][R].sum() * ratio
            if "sim_exp" not in v:
                v["sim_exp"] = refit_spectra(srs)
            Es = [e[R].sum() * ratio for e in v["sim_exp"]]
            lo, med, hi = np.quantile(Es, [0.025, 0.5, 0.975])
            rrows.append((pair, bm, "E", f"peak{u['peak']}", len(srs), f"{e_true:.4g}", f"{lo:.4g}",
                          f"{med:.4g}", f"{hi:.4g}", lo <= e_true <= hi))
        cands = [E_by[m] for m in adm[folded]] + [x for x in (res["E_boot"], res["E_het"], res["E_grid2"])
                                                  if x == x]
        res["E_robust"] = max(cands)
        tier = tier_of(res["E_robust"])
        if res["E_robust"] <= 1 and O_genome > max(5, 3 * res["E_point"]):
            # section 6.3: the empirical count becomes the headline and the tier is read from it
            tier = f"{tier_of(O_genome)} (model fails at the extremes: tier read from O = {O_genome:.0f})"
        if len({tier_of(E_by[m]) for m in adm[folded]}) > 1:
            tier += " [model-dependent]"
        if (pair, bm) in unreliable:
            tier += " [numerically unreliable]"
        res["tier"] = tier
        return res

    for u in units.get(pair, []):
        f_peak = pk["f_blocks"].get((u["chrom"], int(u["pos"]) // 1_000_000), float("nan"))
        res_by = {f: analyse(u, f) for f in (False, True) if VAR[f] is not None}
        # headline per unit (Addenda 1a, 1b): unfolded if the pair has an admissible unfolded model and the
        # index SNP is polarized; otherwise folded (Peak 2 and Peak 6 have no derived-allele delta)
        head = False if (adm[False] and res_by.get(False) is not None) else True
        for folded, res in res_by.items():
            if res is None:
                continue
            L_scan = VAR[folded]["L_scan"]
            ratio = L_scan / VAR[folded]["L_mask"]
            Er, Ep = res["E_robust"], res["E_point"]
            tb = info[res["best"]]["tau_BP"] if res["best"] else float("nan")
            tt = info[res["best"]]["tau_TP"] if res["best"] else float("nan")
            control = pair.startswith("C") or u["own"] != "True"
            urows.append((pair, u["peak"], u["own"], control, "folded" if folded else "unfolded",
                          folded == head, u["f_bp"], u["f_tp"], f"{res['thr']:.4f}",
                          f"{res['O_mask']:.2f}", f"{res['O_genome']:.2f}", res["best"] or "none",
                          f"{Ep / ratio:.4g}", f"{Ep:.4g}",
                          ";".join(f"{m}={e:.4g}" for m, e in res["E_by"].items()),
                          f"{res['E_het']:.4g}", f"{res['E_het_wide']:.4g}", f"{res['E_boot']:.4g}",
                          f"{res['E_grid2']:.4g}", f"{Er:.4g}", res["tier"],
                          # section 4: K_needed = L_scan / E_genome (best model); the E_robust version too
                          f"{L_scan / Ep:.4g}" if Ep == Ep and Ep > 0 else "NA",
                          f"{L_scan / Er:.4g}" if Er == Er and Er > 0 else "NA",
                          f"{tb:.5f}", f"{tt:.5f}", f"{f_peak:.3f}", f"{res['f_star']:.3f}"))

    # tail figure (section 10): genome-wide count of sites at least as extreme as delta, observed outside every
    # GWAS peak and expected under each admissible model, on the projected grid from 0.3 to 1
    for folded in (False, True):
        if VAR[folded] is None:
            continue
        V = VAR[folded]
        ratio = V["L_scan"] / V["L_mask"]
        O_np = sum(np.asarray(arr) for k, arr in V["region_spec"].items() if not k.startswith("peak"))
        for thr in sorted({round(x, 6) for x in np.abs(D).ravel() if x >= 0.3 - EPS}):
            R = sel(folded, thr)
            crows.append((pair, "folded" if folded else "unfolded", f"{thr:.4f}", "observed_nonpeak",
                          f"{O_np[R].sum():.2f}"))
            for m in adm[folded]:
                crows.append((pair, "folded" if folded else "unfolded", f"{thr:.4f}", m,
                              f"{info[m]['exp'][R].sum() * ratio:.4g}"))

    # section 6.4 drift-exposure inputs per pair: tau of the headline analysis's best model, and the
    # model-relative asymmetry of near-fixed differences from the best unfolded model
    for folded in (False, True):
        if not adm[folded]:
            continue
        bm = max(adm[folded], key=lambda m: float(info[m]["r"]["LL"]))
        ub = max(adm[False], key=lambda m: float(info[m]["r"]["LL"])) if adm[False] else None
        asym = info[ub]["asym"] if ub else ["NA"] * 4
        drows.append((pair, "folded" if folded else "unfolded", folded == (not adm[False]), bm,
                      f"{info[bm]['tau_BP']:.5f}", f"{info[bm]['tau_TP']:.5f}",
                      info[bm]["tau_BP"] > info[bm]["tau_TP"], ub or "none", *asym))
    return mrows, lrows, urows, drows, rrows, crows


def do_pair(pair):
    # a SystemExit in a pool worker kills the worker without reporting back; raise it as an error instead
    try:
        return _do_pair(pair)
    except SystemExit as e:
        raise RuntimeError(f"{pair}: {e}") from None


def main():
    global a, best, boot, sims, unreliable, units
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
    boot, sims = collections.defaultdict(list), collections.defaultdict(list)
    for kind, store in (("boot", boot), ("sim", sims)):
        for path in glob.glob(os.path.join(OUT, f"fits/*{a.suffix}.M*.{kind}*.seed*.tsv")):
            with open(path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    # only this run's files: the file is <pair><suffix>.<model>..., so an empty suffix must not
                    # pick up a sensitivity fit such as A1w_ex5.M3... whose rows also say pair A1w
                    if os.path.basename(path).split(".")[0] == r["pair"] + a.suffix:
                        store[(r["pair"], r["model"])].append(r)
    # section 9.3: corner cells that move > 5% when the grid doubles make E "numerically unreliable" (10_checks.py)
    unreliable = set()
    chk = os.path.join(OUT, f"checks/model_checks{a.suffix}.tsv")
    if os.path.exists(chk):
        with open(chk) as fh:
            unreliable = {(r["pair"], r["model"]) for r in csv.DictReader(fh, delimiter="\t")
                          if r["numerically_unreliable"] == "True"}
    else:
        print(f"WARNING no {chk}: section 9.3 numerics check not applied", flush=True)
    units = collections.defaultdict(list)
    with open(os.path.join(OUT, f"tail/index_units{a.suffix}.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            units[r["pair_id"]].append(r)

    # pairs are independent; one worker per pair when SLURM gives several CPUs (output order is unchanged)
    pairs = a.pairs.split(",")
    nproc = min(len(pairs), int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    if nproc > 1:
        with multiprocessing.get_context("fork").Pool(nproc) as pool:
            parts = pool.map(do_pair, pairs, chunksize=1)
    else:
        parts = [do_pair(p) for p in pairs]
    mrows, lrows, urows, drows, rrows, crows = ([r for part in parts for r in part[i]] for i in range(6))

    os.makedirs(os.path.join(OUT, "tail"), exist_ok=True)
    with open(os.path.join(OUT, f"tail/models{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tanalysis\tLL\tconverged\tp_misid\tgate_threshold\tgate_O\tgate_E\tadequate\t"
                "admissible\ttheta\tT_2Nanc\ttau_BP\ttau_TP\tF_BP\tF_TP\tN_anc_mu4e-9\tN_anc_mu1e-8\t"
                "obs_BPder95\texp_BPder95\tobs_TPder95\texp_TPder95\tn_boot\tT_2Nanc_95ci\ttau_BP_95ci\ttau_TP_95ci\t"
                "params\n")
        for r in mrows:
            o.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/ladder{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tthreshold\tdirection\tO\tE\tO_over_E\n")
        for r in lrows:
            o.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/unit_results{a.suffix}.tsv"), "w") as o:
        o.write("pair\tpeak\town\tcontrol\tanalysis\theadline\tf_bp\tf_tp\ttail_threshold\tO_mask\tO_genome\t"
                "best_model\tE_mask_best\tE_genome_best\tE_by_model\tE_hetNe\tE_hetNe_wide_f\tE_boot97.5\tE_grid2\tE_robust\t"
                "tier\tK_needed\tK_needed_robust\ttau_BP\ttau_TP\tf_peak\tf_star\n")
        for r in urows:
            o.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/drift_exposure{a.suffix}.tsv"), "w") as o:
        o.write("pair\tanalysis\theadline\tbest_model\ttau_BP\ttau_TP\ttau_BP_gt_TP\tasym_model\t"
                "obs_BPder95\texp_BPder95\tobs_TPder95\texp_TPder95\n")
        for r in drows:
            o.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/recovery{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tquantity\tunit\tn_sims\ttrue\tq025\tmedian\tq975\tcovered\n")
        for r in rrows:
            o.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(OUT, f"tail/curve{a.suffix}.tsv"), "w") as o:
        o.write("pair\tanalysis\tdelta\tseries\tcount_genome\n")
        for r in crows:
            o.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote tail/models{a.suffix}.tsv, ladder, unit_results ({len(urows)} rows), drift_exposure, "
          f"recovery ({len(rrows)} rows), curve")


if __name__ == "__main__":
    main()
