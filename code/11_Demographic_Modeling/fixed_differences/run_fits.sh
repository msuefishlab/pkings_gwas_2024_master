#!/bin/bash
## run_fits.sh: DESIGN.md section 5. Stage 1 fits every model to the four headline pairs and the two
## TP-TP controls, 20 perturbed starts each (4 tasks x 5 starts), then collects them with Godambe
## uncertainties. Stage 2 (run with "boot" after stage 1 has been collected) refits 100 block-bootstrap
## spectra under M3 and under the highest-likelihood converged model, each started from its best fit.
## Deviation 2 adds refinement rounds between the two stages: 20 more starts per pair x model, perturbed from
## the best start so far (round 2: fold 0.5, seeds 201-204; round 3: fold 0.25, seeds 301-304). "collect"
## then pools every round (08_collect_fits.py --godambe).
##   bash code/11_Demographic_Modeling/fixed_differences/run_fits.sh fit [dependency_jobid]
##   bash code/11_Demographic_Modeling/fixed_differences/run_fits.sh refine 2|3
##   bash code/11_Demographic_Modeling/fixed_differences/run_fits.sh nested 1|2
##   bash code/11_Demographic_Modeling/fixed_differences/run_fits.sh collect [dependency_jobid]
##   bash code/11_Demographic_Modeling/fixed_differences/run_fits.sh boot
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
cd "${root}"
source pkings_gwas.env
dir=code/11_Demographic_Modeling/fixed_differences
logs=output_data/slurm_logs/11_Demographic_Modeling
fd_out=${FD_OUT:-output_data/11_Demographic_Modeling/fixed_differences}
export FD_OUT=${fd_out}
fits=${fd_out}/fits
mkdir -p "${fits}"
PAIRS="A1w A2w A3w A4w C1w C2w"
dx() { echo "singularity exec --bind /mnt/gs21 ${dadi_image} python3 $*"; }

if [[ "${1:-fit}" == "fit" ]]; then
  tasks=${fits}/tasks_fit.txt; : > "${tasks}"
  for p in ${PAIRS}; do for m in ${MODELS:-M0 M1 M2 M3 M0f M1f M2f M3f}; do for s in 1 2 3 4; do
    echo "--pair ${p} --model ${m} --starts 5 --seed ${s} --maxtime 2400" >> "${tasks}"
  done; done; done
  dep=${2:+--dependency=afterok:$2}
  arr=$(sbatch --parsable ${dep} --export=root=${root},FD_OUT=${fd_out},TASKS=${tasks} --array=1-$(wc -l < "${tasks}")%100 ${dir}/07_fit_models.sb)
  col=$(sbatch --parsable --dependency=afterany:${arr} --time=04:00:00 --mem=16G --job-name=fd_collect -o ${logs}/collect_%j.log \
    --wrap "cd ${root} && export OMP_NUM_THREADS=1 && $(dx ${dir}/08_collect_fits.py --godambe)")
  echo "fits ${arr} | collect ${col}"
elif [[ "$1" == "refine" ]]; then
  round=$2
  case ${round} in 2) fold=0.5;; 3) fold=0.25;; *) echo "round 2 or 3"; exit 1;; esac
  tasks=${fits}/tasks_refine${round}.txt; : > "${tasks}"
  # best start so far per pair x model, from the per-start files (every earlier round)
  python3 - "${fits}" ${PAIRS} > "${fits}/refine${round}_starts.txt" <<'PY'
import csv, glob, os, sys
fits, pairs = sys.argv[1], sys.argv[2:]
fixed = {"pair", "model", "target", "seed", "start", "LL", "theta", "seconds"}
for pair in pairs:
    for m in ("M0", "M1", "M2", "M3", "M0f", "M1f", "M2f", "M3f"):
        rows = [r for f in glob.glob(os.path.join(fits, f"{pair}.{m}.data.seed*.tsv"))
                for r in csv.DictReader(open(f), delimiter="\t")]
        if not rows:
            continue
        b = max(rows, key=lambda r: float(r["LL"]))
        print(pair, m, ",".join(b[k] for k in b if k not in fixed))
PY
  while read -r p m p0; do for s in 1 2 3 4; do
    echo "--pair ${p} --model ${m} --starts 5 --seed $((round * 100 + s)) --p0 ${p0} --fold ${fold} --maxtime 2400" >> "${tasks}"
  done; done < "${fits}/refine${round}_starts.txt"
  arr=$(sbatch --parsable --export=root=${root},FD_OUT=${fd_out},TASKS=${tasks} --array=1-$(wc -l < "${tasks}")%150 ${dir}/07_fit_models.sb)
  echo "refine round ${round}: ${arr} ($(wc -l < "${tasks}") tasks)"
elif [[ "$1" == "nested" ]]; then
  # Deviation 2, change 4: a model whose best LL is below its nested sub-model's (M0 < M1 < M2, same folding)
  # gets 20 more starts from the sub-model's best fit, embedded and perturbed by 2^0.25. Pass 1 after round 2
  # (seeds 351-354), pass 2 after round 3 (seeds 361-364).
  pass=$2
  case ${pass} in 1) base=350;; 2) base=360;; *) echo "pass 1 or 2"; exit 1;; esac
  tasks=${fits}/tasks_nested${pass}.txt; : > "${tasks}"
  python3 - "${fits}" ${PAIRS} > "${fits}/nested${pass}_starts.txt" <<'PY'
import csv, glob, os, sys
fits, pairs = sys.argv[1], sys.argv[2:]
fixed = {"pair", "model", "target", "seed", "start", "LL", "theta", "seconds"}
def best(pair, m):
    rows = [r for f in glob.glob(os.path.join(fits, f"{pair}.{m}.data.seed*.tsv"))
            for r in csv.DictReader(open(f), delimiter="\t")]
    return max(rows, key=lambda r: float(r["LL"])) if rows else None
for pair in pairs:
    for f in ("", "f"):
        for sub, sup in (("M0", "M1"), ("M1", "M2")):
            b, B = best(pair, sub + f), best(pair, sup + f)
            if b is None or B is None or float(B["LL"]) >= float(b["LL"]):
                continue
            p = {k: v for k, v in b.items() if k not in fixed}
            if sup == "M2":  # symmetric m -> both directions
                m = p.pop("m")
                p = {**{k: p[k] for k in ("nu_BP", "nu_TP", "T")}, "m_BP_from_TP": m, "m_TP_from_BP": m,
                     **({"p_misid": p["p_misid"]} if "p_misid" in p else {})}
            print(pair, sup + f, ",".join(p.values()), f"{float(b['LL']) - float(B['LL']):.1f}", file=sys.stdout)
PY
  while read -r p m p0 short; do for s in 1 2 3 4; do
    echo "--pair ${p} --model ${m} --starts 5 --seed $((base + s)) --p0 ${p0} --fold 0.25 --maxtime 2400" >> "${tasks}"
  done; done < "${fits}/nested${pass}_starts.txt"
  cat "${fits}/nested${pass}_starts.txt"
  if [[ -s "${tasks}" ]]; then
    arr=$(sbatch --parsable --export=root=${root},FD_OUT=${fd_out},TASKS=${tasks} --array=1-$(wc -l < "${tasks}")%150 ${dir}/07_fit_models.sb)
    echo "nested pass ${pass}: ${arr} ($(wc -l < "${tasks}") tasks)"
  else
    echo "nested pass ${pass}: no violations"
  fi
elif [[ "$1" == "collect" ]]; then
  dep=${2:+--dependency=afterany:$2}
  col=$(sbatch --parsable ${dep} --time=06:00:00 --mem=16G --job-name=fd_collect -o ${logs}/collect_%j.log \
    --wrap "cd ${root} && export OMP_NUM_THREADS=1 && $(dx ${dir}/08_collect_fits.py --godambe)")
  echo "collect ${col}"
else
  tasks=${fits}/tasks_boot.txt; : > "${tasks}"
  # M3 and the best converged model per pair, from best_fits.tsv
  python3 - "${fits}/best_fits.tsv" > "${fits}/boot_targets.txt" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
by = {}
for r in rows:
    by.setdefault(r["pair"], []).append(r)
for pair, rs_all in by.items():
    # unfolded and folded analyses separately: the best converged model of each, plus M3 / M3f
    for folded in (False, True):
        rs = [r for r in rs_all if r["model"].endswith("f") == folded]
        if not rs:
            continue
        conv = [r for r in rs if r["converged"] == "True"] or rs
        top = max(conv, key=lambda r: float(r["LL"]))
        m3 = "M3f" if folded else "M3"
        targets = {top["model"]: top}
        targets.update({r["model"]: r for r in rs if r["model"] == m3})
        for r in targets.values():
            p0 = ",".join(kv.split("=")[1] for kv in r["params"].split(";"))
            print(pair, r["model"], p0)
PY
  while read -r p m p0; do
    # 100 bootstrap spectra per target, five per task (keeps the array under the 1,000-job queue limit)
    for k in $(seq 0 5 95); do
      echo "--pair ${p} --model ${m} --starts 1 --seed $((k + 1000)) --boot $(seq -s, ${k} $((k + 4))) --p0 ${p0} --maxtime 1500" >> "${tasks}"
    done
  done < "${fits}/boot_targets.txt"
  arr=$(sbatch --parsable --export=root=${root},FD_OUT=${fd_out},TASKS=${tasks} --array=1-$(wc -l < "${tasks}")%150 ${dir}/07_fit_models.sb)
  echo "bootstrap refits ${arr} ($(wc -l < "${tasks}") tasks)"
fi
