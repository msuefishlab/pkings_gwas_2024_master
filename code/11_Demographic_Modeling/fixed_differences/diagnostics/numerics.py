#!/usr/bin/env python3
"""Numerics check at the round-1 best fits: negative model cells, roughness and cost for
linear vs log-space extrapolation over grid sizes. Args: pair models(comma)."""
import csv, importlib.util, os, pickle, sys, time, warnings
import numpy as np
import dadi
OUT = os.environ["FD_OUT"]
spec = importlib.util.spec_from_file_location("fitmod", "code/11_Demographic_Modeling/fixed_differences/07_fit_models.py")
fitmod = importlib.util.module_from_spec(spec); spec.loader.exec_module(fitmod)
best = {(r["pair"], r["model"]): r for r in csv.DictReader(open(f"{OUT}/fits/best_fits.tsv"), delimiter="\t")}
pair = sys.argv[1]
pk = pickle.load(open(f"{OUT}/sfs/{pair}.pkl", "rb"))
for model in sys.argv[2].split(","):
    data = fitmod.data_of(pk, model)
    ns = data.sample_sizes
    n = max(ns)
    _, _, pnames, _, _, _ = fitmod.model_setup(model, ns)
    base = fitmod.MODELS[model.rstrip("f")][0]
    f = base if fitmod.is_folded(model) else dadi.Numerics.make_anc_state_misid_func(base)
    p = [float(kv.split("=")[1]) for kv in best[(pair, model)]["params"].split(";")]
    iT = pnames.index("T")
    for ext_name, maker in (("lin", dadi.Numerics.make_extrap_func), ("log", dadi.Numerics.make_extrap_log_func)):
        fex = maker(f)
        for k in (1, 2, 3, 4):
            grid = [k * n + 10, k * n + 20, k * n + 30] if k == 1 else [k * n, k * n + 10 * k, k * n + 20 * k]
            t0 = time.time()
            m = fex(p, ns, grid)
            mf = m.fold() if fitmod.is_folded(model) else m
            arr = np.ma.filled(mf, np.nan)
            dat = np.ma.filled(data, 0.0)
            neg = (arr <= 0) & (dat > 0) & ~np.ma.getmaskarray(data)
            # tail corner cells: delta >= 0.9 (BP-derived), how many of them are negative
            n1, n2 = ns
            I, J = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
            corner = (I / n1 - J / n2) >= 0.9
            lls = []
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for s in range(-5, 6):
                    q = list(p); q[iT] = p[iT] * (1 + 0.002 * s)
                    lls.append(dadi.Inference.ll_multinom(fex(q, ns, grid), data))
            lls = np.array(lls)
            c = np.polyfit(np.arange(-5, 6), lls, 2)
            rough = np.std(lls - np.polyval(c, np.arange(-5, 6)))
            print(f"{pair} {model} {ext_name} grid {grid}: neg cells {neg.sum()} (data {dat[neg].sum():.0f}; "
                  f"corner {int((neg & corner).sum())}) LL {lls[5]:.1f} rough {rough:.2f} "
                  f"({(time.time() - t0) / 12:.2f} s/eval)", flush=True)
