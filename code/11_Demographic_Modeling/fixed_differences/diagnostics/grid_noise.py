#!/usr/bin/env python3
"""Diagnose non-convergence: is the likelihood surface noisy because the grid is too coarse?
For a few best fits: LL at three grid sizes, and LL along a 1% line in T at each grid."""
import csv, importlib.util, os, pickle, sys, time
import numpy as np
import dadi
OUT = os.environ["FD_OUT"]
spec = importlib.util.spec_from_file_location("fitmod", "code/11_Demographic_Modeling/fixed_differences/07_fit_models.py")
fitmod = importlib.util.module_from_spec(spec); spec.loader.exec_module(fitmod)
best = {(r["pair"], r["model"]): r for r in csv.DictReader(open(f"{OUT}/fits/best_fits.tsv"), delimiter="\t")}
for pair, model in [(p, m) for p in sys.argv[1].split(",") for m in sys.argv[2].split(",")]:
    pk = pickle.load(open(f"{OUT}/sfs/{pair}.pkl", "rb"))
    data = fitmod.data_of(pk, model)
    n = max(data.sample_sizes)
    _, _, pnames, _, _, _ = fitmod.model_setup(model, data.sample_sizes)
    base = fitmod.MODELS[model.rstrip("f")][0]
    f = base if fitmod.is_folded(model) else dadi.Numerics.make_anc_state_misid_func(base)
    fex = dadi.Numerics.make_extrap_func(f)
    p = [float(kv.split("=")[1]) for kv in best[(pair, model)]["params"].split(";")]
    iT = pnames.index("T")
    for grid in ([n + 10, n + 20, n + 30], [2 * n, 2 * n + 10, 2 * n + 20], [3 * n, 3 * n + 15, 3 * n + 30]):
        t0 = time.time()
        lls = []
        for k in range(-5, 6):
            q = list(p); q[iT] = p[iT] * (1 + 0.002 * k)
            lls.append(dadi.Inference.ll_multinom(fex(q, data.sample_sizes, grid), data))
        lls = np.array(lls)
        # roughness: residual from a quadratic fit along the line
        c = np.polyfit(np.arange(-5, 6), lls, 2)
        rough = np.std(lls - np.polyval(c, np.arange(-5, 6)))
        print(f"{pair} {model} grid {grid}: LL@best {lls[5]:.2f}  range {lls.max()-lls.min():.2f}  "
              f"roughness(sd resid quad) {rough:.3f}  ({time.time()-t0:.0f} s for 11 evals)", flush=True)
