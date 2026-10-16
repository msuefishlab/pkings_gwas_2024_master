#!/usr/bin/env python3
"""13_report_tables.py
Markdown tables for RESULTS.md from the 09_tail.py outputs (DESIGN.md section 10). Standard library only.

  python3 code/11_Demographic_Modeling/fixed_differences/13_report_tables.py <tail_dir> > tables.md

Sections: own units (headline analysis), the sentence of section 10 for each own unit, off-target and control
units by tier, models, drift exposure (section 6.4), and parameter recovery (section 9.3).
"""
import collections
import csv
import math
import os
import sys

POP = {"A1w": ("BAM", "APA + BEN"), "A2w": ("BIK", "DOV + DOG"), "A3w": ("BAVA", "DOV + DOG"),
       "A4w": ("BIR", "DOV + DOG"), "C1w": ("APA", "BEN"), "C2w": ("DOV", "DOG")}


def read(tail, name):
    with open(os.path.join(tail, name)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def fmt(x, nd=3):
    v = num(x)
    if v != v:
        return "NA"
    if v == 0:
        return "0"
    if abs(v) >= 1e5 or abs(v) < 1e-3:
        return f"{v:.{nd - 1}e}"
    if abs(v) >= 100:
        return f"{v:,.0f}"
    return f"{v:.{nd}g}"


def share(e):
    # section 2: K_needed as a fraction of all callable sites, 1 / E, at the drift-generous E_robust
    v = num(e)
    if v != v:
        return "NA"
    if v <= 1:
        return "> 100%"
    x = 100 / v
    return f"{x:.2g}%" if x < 1 else f"{x:.0f}%"


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def main():
    tail = sys.argv[1]
    units = read(tail, "unit_results.tsv")
    models = {(r["pair"], r["model"]): r for r in read(tail, "models.tsv")}
    drift = read(tail, "drift_exposure.tsv")
    recov = read(tail, "recovery.tsv") if os.path.exists(os.path.join(tail, "recovery.tsv")) else []

    own = [u for u in units if u["own"] == "True" and u["headline"] == "True" and not u["pair"].startswith("C")]
    print("## Own units (headline analysis)\n")
    rows = []
    for u in own:
        ci = (f"{fmt(u['E_boot2.5'])}–{fmt(u['E_boot97.5'])}" if num(u.get("E_boot2.5")) == num(u.get("E_boot2.5"))
              else "NA")
        rows.append((f"{POP[u['pair']][0]} vs {POP[u['pair']][1]}", f"P{u['peak']}", u["analysis"],
                     fmt(u["tail_threshold"]), fmt(u["O_genome"]), u["best_model"], fmt(u["E_genome_best"]), ci,
                     fmt(u["E_hetNe"]), fmt(u["E_grid2"]), fmt(u["E_robust"]), u["tier"], fmt(u["K_needed"]),
                     share(u["E_robust"])))
    print(table(["Pair", "Peak", "Analysis", "Δ", "O (non-peak)", "Best model", "E", "E 95% CI", "E_hetNe",
                 "E_grid2", "E_robust", "Tier", "K_needed", "BP-capable share needed"], rows))
    print("\nK_needed = L_scan / E (best model). The last column is 1 / E_robust: the share of all callable sites "
          "that would each have to be able to make a fish BP for drift to be expected to produce one site this "
          "differentiated (section 2), at the most drift-generous E.")

    print("\n## Section 10 sentences\n")
    for u in own:
        bp, tp = POP[u["pair"]]
        bm = u["best_model"]
        if bm == "none":
            print(f"- **{bp}, P{u['peak']}** ({u['analysis']}): no admissible model, so no E. We observe "
                  f"{fmt(u['O_genome'])} sites outside the GWAS peaks at least as differentiated (|Δ| ≥ "
                  f"{fmt(u['tail_threshold'])}). Verdict: {u['tier']}.")
            continue
        m = models[(u["pair"], bm)]
        F, t2n, tau = num(m["F_BP"]), num(m["T_2Nanc"]), num(m["tau_BP"])
        # A new neutral mutation that fixes takes on average 4N generations, tau = 2 in units of 2N. The clock is
        # the BP branch's own drift (tau_BP), not T/2N_anc: the BP population is smaller than the ancestor in every
        # headline fit, so tau_BP is the larger, drift-generous value.
        frac = tau / 2 * 100
        if frac < 25:
            new_ok = "so a neutral fixed difference must come from variation that already existed"
        elif frac < 100:
            new_ok = ("so new mutations could only rarely have drifted to fixation; most neutral fixed differences "
                      "come from variation that already existed")
        else:
            new_ok = "enough for new neutral mutations to have drifted to fixation"
        ci = (f" (bootstrap 95% interval {fmt(u['E_boot2.5'])}–{fmt(u['E_boot97.5'])})"
              if num(u.get("E_boot2.5")) == num(u.get("E_boot2.5")) else "")
        print(f"- **{bp}, P{u['peak']}** ({u['analysis']}, {bm}): Since {bp} split from {tp}, drift has had time "
              f"to move a typical 50% allele by about ±{math.sqrt(F * 0.25):.2f}. The elapsed time is "
              f"{t2n:.3g} × 2N_anc generations; in {bp}'s own drift units it is τ_BP = {tau:.3g}, {frac:.0f}% of "
              f"the time a new neutral mutation needs, on average, to drift to fixation, {new_ok}. Under the best-fitting history, drift would create {fmt(u['E_genome_best'])} sites{ci} "
              f"in the genome as strongly BP-biased as this SNP (most drift-generous value, E_robust: "
              f"{fmt(u['E_robust'])}). We observe {fmt(u['O_genome'])} outside the GWAS peaks. "
              f"Verdict: {u['tier']}.")

    print("\n## Off-target and control units (headline analysis)\n")
    cnt = collections.Counter()
    for u in units:
        if u["headline"] != "True" or u in own:
            continue
        kind = "control pair" if u["pair"].startswith("C") else "off-target"
        cnt[(kind, u["tier"].split(" [")[0].split(" (")[0])] += 1
    print(table(["Unit kind", "Tier", "Units"], [(k, t, n) for (k, t), n in sorted(cnt.items())]))

    print("\n## Models\n")
    rows = []
    for (pair, model), m in models.items():
        rows.append((pair, model, m["converged"], m["admissible"], fmt(m["p_misid"]), fmt(m["T_2Nanc"]),
                     m["T_2Nanc_95ci"], fmt(m["tau_BP"]), fmt(m["tau_TP"]), m["tau_BP_95ci"], m["tau_TP_95ci"],
                     m["n_boot"], m.get("at_bound", "NA")))
    print(table(["Pair", "Model", "Converged", "Admissible", "p_misid", "T/2N_anc", "95% CI", "τ_BP", "τ_TP",
                 "τ_BP 95% CI", "τ_TP 95% CI", "Bootstrap refits", "At a bound"], rows))

    print("\n## Drift exposure (section 6.4)\n")
    rows = [(f"{POP[d['pair']][0]} vs {POP[d['pair']][1]}", d["analysis"], d["headline"], d["best_model"],
             fmt(d["tau_BP"]), fmt(d["tau_TP"]), d["tau_BP_gt_TP"], fmt(d["obs_BPder95"]), fmt(d["exp_BPder95"]),
             fmt(d["obs_TPder95"]), fmt(d["exp_TPder95"])) for d in drift]
    print(table(["Pair", "Analysis", "Headline", "Best model", "τ_BP", "τ_TP", "τ_BP > τ_TP",
                 "BP-derived ≥0.95 obs", "exp", "TP-derived ≥0.95 obs", "exp"], rows))
    head = [d for d in drift if d["headline"] == "True" and not d["pair"].startswith("C")]
    k = sum(d["tau_BP_gt_TP"] == "True" for d in head)
    print(f"\nτ_BP > τ_TP in {k} of {len(head)} headline pairs.")

    print("\n## Parameter recovery (section 9.3)\n")
    if not recov:
        print("No Poisson refits yet.")
    else:
        rows = [(r["pair"], r["model"], r["quantity"], r["unit"], r["n_sims"], fmt(r["true"]),
                 f"{fmt(r['q025'])}–{fmt(r['q975'])}", r["covered"]) for r in recov]
        print(table(["Pair", "Model", "Quantity", "Unit", "Replicates", "True", "95% range", "Covered"], rows))


if __name__ == "__main__":
    main()
