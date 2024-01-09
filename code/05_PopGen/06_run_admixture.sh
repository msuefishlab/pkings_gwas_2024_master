#!/bin/bash

## 06_run_admixture.sh
## Jan 2024
## JRG
# This script runs ADMIXTURE on the GDS file genereated for the phylogeny

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

MYNAME=$(whoami)

indir=${root}/output_data/04_Phylogeny/
outdir=${root}/output_data/05_PopGen/

mkdir -p ${outdir}
mkdir -p ${tmpoutdir}

OUTROOT=${snps_only_vcf%".vcf.gz"}


MYNAME=$(whoami)

singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} R --slave --vanilla --file="/project_root/code/05_PopGen/convert_gds_to_plink.R" --args -i "/indir/${OUTROOT}.renamed.maf10.miss10.dp5.gds" -o "/outdir/${OUTROOT}.renamed.maf10.miss10.dp5.filtered"

cat ${outdir}/{OUTROOT}.renamed.maf10.miss10.dp5.filtered.fam | cut -f2 | awk 'BEGIN { FS = "_" } ; { print $1"_"$2,"unk",$1 }' > ${outdir}/indfile.txt

numk=$(cat ${outdir}/indfile.txt | cut -d" " -f3 | uniq | wc -l)

for K in $(seq 1 $numk); do
    CMD="singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir ${gwas_tools_image} /admixture/admixture --cv -j8 /out_dir/${OUTROOT}.renamed.maf10.miss10.dp5.filtered.bed $K | tee /out_dir/{OUTROOT}.renamed.maf10.miss10.dp5.filtered.${K}.out" #Normally you should give --cv as first option to admixture
    echo sbatch --time=04:00:00 -c 8 -J ADMIXTURE_${K} -o ${root}"/output_data/slurm_logs/05_PopGen/popgen_ADMIXTURE_${K}.log" --export=root=${root} --mem 12000 --wrap="$CMD"
    sbatch --time=04:00:00 -c 8 -J ADMIXTURE_${K} -o ${root}"/output_data/slurm_logs/05_PopGen/popgen_ADMIXTURE_${K}.log" --export=root=${root} --mem 12000 --wrap="$CMD"
done
