#!/usr/bin/env python3
"""08_collect_fits.py
Collect the starts from 07_fit_models.py, check convergence, and add Godambe uncertainties and the
adjusted likelihood-ratio tests (DESIGN.md section 5).

Converged: the best three starts lie within 1 log-likelihood unit and their parameters within 5%.
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

OUT = os.path.join(os.getcwd(), "output_data/11_Demographic_Modeling/fixed_differences")
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
    for path in glob.glob(os.path.join(OUT, f"fits/*{a.suffix}.M?.data.seed*.tsv")):
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
        row = {"pair": pair, "model": model, "n_starts": len(rs), "LL": float(r0["LL"]),
               "theta": float(r0["theta"]), "converged": ll_ok and par_ok,
               "LL_gap_top3": float(top[0]["LL"]) - float(top[-1]["LL"]) if len(top) >= 3 else float("nan"),
               "params": {p: float(r0[p]) for p in pnames}}
        if a.godambe:
            with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.pkl"), "rb") as fh:
                pk = pickle.load(fh)
            data = dadi.Spectrum(pk["data"], mask_corners=True)
            func_ex, pts_l, _, _, _, fixed = fitmod.model_setup(model, data.sample_sizes)
            boots = fitmod.bootstrap_spectra(pk, block_mb=a.block_mb)
            free = [i for i, f in enumerate(fixed or [None] * len(pnames)) if f is None]
            popt = [row["params"][p] for p in pnames]

            def f_free(pf, ns, pts):
                full = list(popt)
                for i, v in zip(free, pf):
                    full[i] = v
                return func_ex(full, ns, pts)
            try:
                se = dadi.Godambe.GIM_uncert(f_free, pts_l, boots, [popt[i] for i in free], data, multinom=True)
                row["se"] = {pnames[i]: float(s) for i, s in zip(free, se[:-1])}
                row["se_theta"] = float(se[-1])
            except Exception as e:                       # report, do not hide
                row["se"] = {}
                row["godambe_error"] = str(e)
        rows.append(row)
        best[(pair, model)] = row
        print(pair, model, f"LL {row['LL']:.2f}", "converged" if row["converged"] else "NOT converged",
              " ".join(f"{p}={v:.4g}" for p, v in row["params"].items()))

    with open(os.path.join(OUT, f"fits/best_fits{a.suffix}.tsv"), "w") as o:
        o.write("pair\tmodel\tn_starts\tLL\ttheta\tconverged\tLL_gap_top3\tparams\tse\n")
        for r in rows:
            o.write(f"{r['pair']}\t{r['model']}\t{r['n_starts']}\t{r['LL']:.4f}\t{r['theta']:.6g}\t{r['converged']}\t"
                    f"{r['LL_gap_top3']:.3f}\t"
                    + ";".join(f"{k}={v:.6g}" for k, v in r["params"].items()) + "\t"
                    + ";".join(f"{k}={v:.4g}" for k, v in r.get("se", {}).items()) + "\n")

    # adjusted LRTs (information only): M0 (no migration) nested in M1 and in M2
    if a.godambe:
        with open(os.path.join(OUT, f"fits/lrt{a.suffix}.tsv"), "w") as o:
            o.write("pair\tsimple\tcomplex\tD_raw\tD_adj\tp_adj\n")
            for pair in sorted({p for p, _ in best}):
                for comp, nested in (("M1", [3]), ("M2", [3, 4])):
                    if (pair, "M0") not in best or (pair, comp) not in best:
                        continue
                    with open(os.path.join(OUT, f"sfs/{pair}{a.suffix}.pkl"), "rb") as fh:
                        pk = pickle.load(fh)
                    data = dadi.Spectrum(pk["data"], mask_corners=True)
                    func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(comp, data.sample_sizes)
                    boots = fitmod.bootstrap_spectra(pk, block_mb=a.block_mb)
                    p_c = [best[(pair, comp)]["params"][p] for p in pnames]
                    p0 = list(p_c)
                    for i in nested:
                        p0[i] = 0.0
                    p0[:3] = [best[(pair, "M0")]["params"][p] for p in pnames[:3]]
                    p0[-1] = best[(pair, "M0")]["params"]["p_misid"]
                    d_raw = 2 * (best[(pair, comp)]["LL"] - best[(pair, "M0")]["LL"])
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
                    o.write(f"{pair}\tM0\t{comp}\t{d_raw:.3f}\t{d_adj:.3f}\t{p:.4g}\n")


if __name__ == "__main__":
    main()
