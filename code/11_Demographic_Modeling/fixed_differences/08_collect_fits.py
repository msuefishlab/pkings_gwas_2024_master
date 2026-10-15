#!/usr/bin/env python3
"""08_collect_fits.py
Collect the starts from 07_fit_models.py, check convergence, and add Godambe uncertainties and the
adjusted likelihood-ratio tests (DESIGN.md section 5).

Converged (DESIGN.md Deviation 2): the best three starts, pooled over all rounds, lie within one
Godambe-adjusted log-likelihood unit: (LL_1 - LL_3) / c <= 1, where c is the mean eigenvalue of H^-1 J at the
best start (the composite-likelihood inflation; GIM = H J^-1 H). Needs --godambe; without it the strict rule of
section 5 is used (best three within 1 raw unit, parameters within 5%), which is also always reported.
Godambe (GIM) uses 100 bootstrap spectra over 1-Mb blocks (--block_mb 5 for the sensitivity row).
LRTs are reported for information only (M0 nested in M1 and in M2, migration fixed at 0); no AIC,
because dadi's likelihood is composite.

Writes fits/best_fits<suffix>.tsv (one row per pair x model) and fits/lrt<suffix>.tsv.
"""
import argparse
import collections
import csv
import glob
import importlib.util
import os
import pickle

import numpy as np
import dadi

OUT = os.path.join(os.getcwd(), os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("fitmod", os.path.join(HERE, "07_fit_models.py"))
fitmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitmod)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="")
    ap.add_argument("--godambe", action="store_true")
    ap.add_argument("--block_mb", type=int, default=1)
    a = ap.parse_args()

    runs = collections.defaultdict(list)
    for path in glob.glob(os.path.join(OUT, f"fits/*{a.suffix}.M*.data.seed*.tsv")):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                runs[(r["pair"], r["model"])].append(r)

    rows, best = [], {}
    for (pair, model), rs in sorted(runs.items()):
        rs.sort(key=lambda r: -float(r["LL"]))
        pnames = [k for k in rs[0] if k not in ("pair", "model", "target", "seed", "start", "LL", "theta", "seconds")]
        top = rs[:3]
        p_top = np.array([[float(r[p]) for p in pnames] for r in top])
        ll_ok = len(top) >= 3 and float(top[0]["LL"]) - float(top[-1]["LL"]) <= 1.0
        rel = np.abs(p_top - p_top[0]) / np.maximum(np.abs(p_top[0]), 1e-8)
        # parameters pinned near zero (e.g. migration) count as agreeing if all are below 1e-3
        small = np.all(np.abs(p_top) < 1e-3, axis=0)
        par_ok = bool(np.all((rel[1:] <= 0.05) | small))
        r0 = rs[0]
        gap = float(top[0]["LL"]) - float(top[-1]["LL"]) if len(top) >= 3 else float("nan")
        row = {"pair": pair, "model": model, "n_starts": len(rs), "LL": float(r0["LL"]),
               "theta": float(r0["theta"]), "converged_strict": ll_ok and par_ok, "converged": ll_ok and par_ok,
               "LL_gap_top3": gap, "c_adj": float("nan"), "gap_adj": float("nan"),
               "params": {p: float(r0[p]) for p in pnames},
               "top3": [{p: float(r[p]) for p in pnames} for r in top]}
        if a.godambe:
            with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.pkl"), "rb") as fh:
                pk = pickle.load(fh)
            data = fitmod.data_of(pk, model)
            func_ex, pts_l, _, _, _, fixed = fitmod.model_setup(model, data.sample_sizes)
            boots = fitmod.bootstrap_spectra(pk, block_mb=a.block_mb, folded=fitmod.is_folded(model))
            free = [i for i, f in enumerate(fixed or [None] * len(pnames)) if f is None]
            popt = [row["params"][p] for p in pnames]

            def f_free(pf, ns, pts):
                full = list(popt)
                for i, v in zip(free, pf):
                    full[i] = v
                return func_ex(full, ns, pts)
            try:
                se, gim, hess = dadi.Godambe.GIM_uncert(f_free, pts_l, boots, [popt[i] for i in free], data,
                                                        multinom=True, return_GIM=True)
                row["se"] = {pnames[i]: float(s) for i, s in zip(free, se[:-1])}
                row["se_theta"] = float(se[-1])
                # H^-1 J = GIM^-1 H; its mean eigenvalue is the factor by which linkage inflates the raw
                # composite log-likelihood (Deviation 2). A non-positive value (best start not at an optimum)
                # leaves the model unconverged.
                c = float(np.mean(np.real(np.linalg.eigvals(np.linalg.inv(gim) @ hess))))
                row["c_adj"] = c
                if c > 0 and len(top) >= 3:
                    row["gap_adj"] = gap / c
                    row["converged"] = gap / c <= 1.0
                else:
                    row["converged"] = False
            except Exception as e:                       # report, do not hide
                row["se"] = {}
                row["godambe_error"] = str(e)
                row["converged"] = False
        rows.append(row)
        best[(pair, model)] = row
        print(pair, model, f"LL {row['LL']:.2f} gap {gap:.1f} c {row['c_adj']:.4g} gap_adj {row['gap_adj']:.3g}",
              "converged" if row["converged"] else "NOT converged",
              "(strict yes)" if row["converged_strict"] else "(strict no)",
              " ".join(f"{p}={v:.4g}" for p, v in row["params"].items()), flush=True)

    with open(os.path.join(OUT, f"fits/best_fits{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tn_starts\tLL\ttheta\tconverged\tLL_gap_top3\tparams\tse\t"
                "converged_strict\tc_adj\tgap_adj\ttop3_params\n")
        for r in rows:
            o.write(f"{r['pair']}\t{r['model']}\t{r['n_starts']}\t{r['LL']:.4f}\t{r['theta']:.6g}\t{r['converged']}\t"
                    f"{r['LL_gap_top3']:.3f}\t"
                    + ";".join(f"{k}={v:.6g}" for k, v in r["params"].items()) + "\t"
                    + ";".join(f"{k}={v:.4g}" for k, v in r.get("se", {}).items()) + "\t"
                    + f"{r['converged_strict']}\t{r['c_adj']:.6g}\t{r['gap_adj']:.4g}\t"
                    + "|".join(";".join(f"{k}={v:.8g}" for k, v in t.items()) for t in r["top3"]) + "\n")

    # adjusted LRTs (information only): M0 (no migration) nested in M1 and in M2
    if a.godambe:
        with open(os.path.join(OUT, f"fits/lrt{a.suffix}.tsv"), "w") as o:
            o.write("pair\tsimple\tcomplex\tD_raw\tD_adj\tp_adj\n")
            for pair in sorted({p for p, _ in best}):
                for simple, comp, nested in (("M0", "M1", [3]), ("M0", "M2", [3, 4]),
                                             ("M0f", "M1f", [3]), ("M0f", "M2f", [3, 4])):
                    if (pair, simple) not in best or (pair, comp) not in best:
                        continue
                    with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.pkl"), "rb") as fh:
                        pk = pickle.load(fh)
                    data = fitmod.data_of(pk, comp)
                    func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(comp, data.sample_sizes)
                    boots = fitmod.bootstrap_spectra(pk, block_mb=a.block_mb, folded=fitmod.is_folded(comp))
                    p_c = [best[(pair, comp)]["params"][p] for p in pnames]
                    p0 = list(p_c)
                    for i in nested:
                        p0[i] = 0.0
                    p0[:3] = [best[(pair, simple)]["params"][p] for p in pnames[:3]]
                    if "p_misid" in pnames:
                        p0[-1] = best[(pair, simple)]["params"]["p_misid"]
                    d_raw = 2 * (best[(pair, comp)]["LL"] - best[(pair, simple)]["LL"])
                    try:
                        adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, boots, p0, data, nested_indices=nested, multinom=True)
                        d_adj = adj * d_raw
                        # migration fixed at its boundary (0): chi-square mixtures (Self and Liang 1987)
                        from scipy.stats import chi2
                        if len(nested) == 1:
                            p = 0.5 * chi2.sf(d_adj, 1)
                        else:
                            p = 0.5 * chi2.sf(d_adj, 1) + 0.25 * chi2.sf(d_adj, 2)
                    except Exception as e:
                        d_adj, p = float("nan"), float("nan")
                    o.write(f"{pair}\t{simple}\t{comp}\t{d_raw:.3f}\t{d_adj:.3f}\t{p:.4g}\n")


if __name__ == "__main__":
    main()
