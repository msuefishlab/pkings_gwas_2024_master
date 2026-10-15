#!/usr/bin/env python3
"""07_fit_models.py
Fit one two-population model to one pair's joint spectrum from several perturbed starts (DESIGN.md
section 5). BP is population 1 throughout, so in IM the founder fraction s belongs to BP.

  M0  split, no migration          split_mig with m fixed at 0       (nu1, nu2, T, m=0)
  M1  split, symmetric migration   split_mig                         (nu1, nu2, T, m)
  M2  split, asymmetric migration  split_asym_mig                    (nu1, nu2, T, m12, m21)
  M3  isolation with migration     IM: founder fraction s, growth    (s, nu1, nu2, T, m12, m21)
Every model is wrapped with dadi's ancestral-misidentification parameter (last parameter).

Writes one row per start to fits/<pair><suffix>.<model>.seed<seed>.tsv.
With --boot K, fits bootstrap spectrum K instead of the data, starting from --p0 (for confidence intervals).

Usage (inside dadi_tools.sif, from the repo root):
  python3 .../07_fit_models.py --pair A1w --model M3 --starts 5 --seed 1
"""
import argparse
import os
import pickle
import time

import numpy as np
import dadi
import nlopt

OUT = os.path.join(os.getcwd(), "output_data/11_Demographic_Modeling/fixed_differences")

MODELS = {
    "M0": (dadi.Demographics2D.split_mig, ["nu_BP", "nu_TP", "T", "m"], [None, None, None, 0.0]),
    "M1": (dadi.Demographics2D.split_mig, ["nu_BP", "nu_TP", "T", "m"], None),
    "M2": (dadi.Demographics2D.split_asym_mig, ["nu_BP", "nu_TP", "T", "m_BP_from_TP", "m_TP_from_BP"], None),
    "M3": (dadi.Demographics2D.IM, ["s_BP", "nu_BP", "nu_TP", "T", "m_BP_from_TP", "m_TP_from_BP"], None),
}
BOUNDS = {"nu_BP": (1e-3, 100), "nu_TP": (1e-3, 100), "T": (1e-5, 10), "m": (0, 30),
          "m_BP_from_TP": (0, 30), "m_TP_from_BP": (0, 30), "s_BP": (1e-3, 0.999), "p_misid": (0, 0.2)}
START = {"nu_BP": 1.0, "nu_TP": 1.0, "T": 0.2, "m": 1.0, "m_BP_from_TP": 1.0, "m_TP_from_BP": 1.0,
         "s_BP": 0.3, "p_misid": 0.02}


def model_setup(name, ns):
    base, pnames, fixed = MODELS[name]
    pnames = pnames + ["p_misid"]
    func = dadi.Numerics.make_anc_state_misid_func(base)
    func_ex = dadi.Numerics.make_extrap_func(func)
    n = max(ns)
    pts_l = [n + 10, n + 20, n + 30]
    lower = [BOUNDS[p][0] for p in pnames]
    upper = [BOUNDS[p][1] for p in pnames]
    fixed = (fixed + [None]) if fixed else None
    return func_ex, pts_l, pnames, lower, upper, fixed


def bootstrap_spectra(pk, n_boot=100, seed=20261015, block_mb=1):
    """Block bootstrap over 1-Mb blocks (or unions of block_mb consecutive blocks on one chromosome)."""
    keys, stack = pk["block_keys"], pk["block_spec"]
    groups = {}
    for i, (c, b) in enumerate(keys):
        groups.setdefault((c, b // block_mb), []).append(i)
    units = [np.sum(stack[idx], axis=0) for idx in groups.values()]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(units), len(units))
        out.append(dadi.Spectrum(np.sum([units[i] for i in pick], axis=0), mask_corners=True))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", required=True)
    ap.add_argument("--model", required=True, choices=sorted(MODELS))
    ap.add_argument("--starts", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--boot", default=None, help="refit bootstrap spectra K[,K2,...] (needs --p0)")
    ap.add_argument("--p0", default=None, help="comma-separated start parameters (bootstrap refits)")
    ap.add_argument("--maxtime", type=float, default=2400, help="seconds per start")
    a = ap.parse_args()

    with open(os.path.join(OUT, f"sfs/{a.pair}{a.suffix}.pkl"), "rb") as fh:
        pk = pickle.load(fh)
    targets = [("data", dadi.Spectrum(pk["data"], mask_corners=True))]
    if a.boot is not None:
        boots = bootstrap_spectra(pk)
        targets = [(f"boot{k}", boots[int(k)]) for k in a.boot.split(",")]
    os.makedirs(os.path.join(OUT, "fits"), exist_ok=True)
    np.random.seed(a.seed)
    for tag, data in targets:
        fit_target(a, tag, data)


def fit_target(a, tag, data):
    ns = data.sample_sizes
    func_ex, pts_l, pnames, lower, upper, fixed = model_setup(a.model, ns)
    out_path = os.path.join(OUT, f"fits/{a.pair}{a.suffix}.{a.model}.{tag}.seed{a.seed}.tsv")
    with open(out_path, "w") as o:
        o.write("\t".join(["pair", "model", "target", "seed", "start", "LL", "theta", "seconds"] + pnames) + "\n")
        for k in range(a.starts):
            if a.p0:
                p0 = [float(v) for v in a.p0.split(",")]
                p0 = dadi.Misc.perturb_params(p0, fold=0.5, lower_bound=lower, upper_bound=upper)
            else:
                p0 = dadi.Misc.perturb_params([START[p] for p in pnames], fold=1.5,
                                              lower_bound=lower, upper_bound=upper)
            if fixed:
                p0 = [f if f is not None else v for v, f in zip(p0, fixed)]
            t0 = time.time()
            res = dadi.Inference.opt(p0, data, func_ex, pts_l, lower_bound=lower, upper_bound=upper,
                                     fixed_params=fixed, algorithm=nlopt.LN_BOBYQA, maxtime=a.maxtime)
            popt = res[0] if isinstance(res, tuple) else res
            model = func_ex(popt, ns, pts_l)
            ll = dadi.Inference.ll_multinom(model, data)
            theta = dadi.Inference.optimal_sfs_scaling(model, data)
            o.write("\t".join([a.pair, a.model, tag, str(a.seed), str(k), f"{ll:.4f}", f"{theta:.6g}",
                               f"{time.time() - t0:.0f}"] + [f"{v:.8g}" for v in popt]) + "\n")
            o.flush()
            print(f"{tag} start {k}: LL {ll:.2f} theta {theta:.4g} ({time.time() - t0:.0f} s) "
                  + " ".join(f"{p}={v:.4g}" for p, v in zip(pnames, popt)), flush=True)


if __name__ == "__main__":
    main()
