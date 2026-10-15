#!/bin/bash
## run_to_spectra.sh: chain everything from the finished site table to the gathered spectra.
## The gate (02_check_site_table.py) runs first; it exits non-zero on failure, and every later job
## depends on it with afterok, so nothing downstream runs on a site table that failed.
##   bash code/11_Demographic_Modeling/fixed_differences/run_to_spectra.sh [site_table_jobid]
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
cd "${root}"
source pkings_gwas.env
dir=code/11_Demographic_Modeling/fixed_differences
logs=output_data/slurm_logs/11_Demographic_Modeling
n=$(( $(wc -l < input_data/11_Demographic_Modeling/chunks.tsv) - 2 ))
dep=${1:+--dependency=afterok:$1}
dx() { echo "singularity exec --bind /mnt/gs21 ${dadi_image} python3 $*"; }
fd_out=${FD_OUT:-output_data/11_Demographic_Modeling/fixed_differences}
export FD_OUT=${fd_out}
mkdir -p ${fd_out}/sfs/chunks

gate=$(sbatch --parsable ${dep} --time=01:00:00 --mem=4G --job-name=fd_gate -o ${logs}/gate_%j.log \
  --wrap "cd ${root} && $(dx ${dir}/02_check_site_table.py)")
anh=$(sbatch --parsable --dependency=afterok:${gate} --export=root=${root},FD_OUT=${fd_out},MODE=an_hist --array=0-${n}%60 ${dir}/03_spectra.sb)
proj=$(sbatch --parsable --dependency=afterok:${anh} --time=00:20:00 --mem=2G --job-name=fd_proj -o ${logs}/proj_%j.log \
  --wrap "cd ${root} && $(dx ${dir}/04_projections.py)")
bins=$(sbatch --parsable --dependency=afterok:${proj} --export=root=${root},FD_OUT=${fd_out},MODE=bins --array=0-${n}%60 ${dir}/03_spectra.sb)
des=$(sbatch --parsable --dependency=afterok:${bins} --time=00:20:00 --mem=4G --job-name=fd_deserts -o ${logs}/deserts_%j.log \
  --wrap "cd ${root} && $(dx ${dir}/05_deserts.py)")
spec=$(sbatch --parsable --dependency=afterok:${des} --export=root=${root},FD_OUT=${fd_out},MODE=spectra --array=0-${n}%60 ${dir}/03_spectra.sb)
gat=$(sbatch --parsable --dependency=afterok:${spec} --time=01:00:00 --mem=16G --job-name=fd_gather -o ${logs}/gather_%j.log \
  --wrap "cd ${root} && $(dx ${dir}/06_gather.py)")
echo "gate ${gate} | an_hist ${anh} | projections ${proj} | bins ${bins} | deserts ${des} | spectra ${spec} | gather ${gat}"
