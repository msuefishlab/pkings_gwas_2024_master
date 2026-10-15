#!/usr/bin/env python3
"""Roughness of the LL along T at the round-1 best fit, by time step (timescale_factor) and grid.
Args: pair model."""
import csv, importlib.util, os, pickle, sys, time, warnings
import numpy as np
import dadi
OUT = os.environ["FD_OUT"]
spec = importlib.util.spec_from_file_location("fitmod", "code/11_Demographic_Modeling/fixed_differences/07_fit_models.py")
fitmod = importlib.util.module_from_spec(spec); spec.loader.exec_module(fitmod)
best = {(r["pair"], r["model"]): r for r in csv.DictReader(open(f"{OUT}/fits/best_fits.tsv"), delimiter="\t")}
pair, model = sys.argv[1], sys.argv[2]
pk = pickle.load(open(f"{OUT}/sfs/{pair}.pkl", "rb"))
data = fitmod.data_of(pk, model)
ns = data.sample_sizes; n = max(ns)
_, _, pnames, _, _, _ = fitmod.model_setup(model, ns)
base = fitmod.MODELS[model.rstrip("f")][0]
f = base if fitmod.is_folded(model) else dadi.Numerics.make_anc_state_misid_func(base)
fex = dadi.Numerics.make_extrap_log_func(f)
p = [float(kv.split("=")[1]) for kv in best[(pair, model)]["params"].split(";")]
iT = pnames.index("T")
for tsf in (1e-3, 2.5e-4, 1e-4):
    dadi.Integration.timescale_factor = tsf
    for k in (1, 2):
        grid = [n + 10, n + 20, n + 30] if k == 1 else [2 * n, 2 * n + 20, 2 * n + 40]
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lls = np.array([dadi.Inference.ll_multinom(fex([*p[:iT], p[iT] * (1 + 0.002 * s), *p[iT + 1:]], ns, grid), data)
                            for s in range(-5, 6)])
        c = np.polyfit(np.arange(-5, 6), lls, 2)
        rough = np.std(lls - np.polyval(c, np.arange(-5, 6)))
        print(f"{pair} {model} tsf {tsf:g} grid {grid}: LL {lls[5]:.1f} rough {rough:.2f} ({(time.time()-t0)/11:.2f} s/eval)", flush=True)
