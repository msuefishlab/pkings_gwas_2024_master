#!/bin/bash

## 03_run_admixture.sh
## Jan 2024
## JRG
# This script runs ADMIXTURE on the pruned SNP set for a range of K values
# Usage: bash 02_run_admixture.sh OUTNAME NUMK

root="$(git rev-parse --show-toplevel)"

OUTNAME=$1

NUMK=$2  # Number of populations to test

source $root/pkings_gwas.env

MYNAME=$(whoami)

indir=${root}/output_data/13_Structure/
keepdir=${root}/output_data/13_Structure/$OUTNAME/
outdir=${root}/output_data/13_Structure/$OUTNAME/

mkdir -p ${outdir}

OUTROOT=${snps_only_vcf%".vcf.gz"}

singularity exec --bind $root:/project_root --bind $keepdir:/keep_dir --bind $outdir:/out_dir --bind $indir:/in_dir ${gwas_tools_image} /plink/plink \
--bfile /in_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture  \
--allow-extra-chr \
--allow-no-sex \
--make-bed \
--keep /keep_dir/${OUTNAME}.keep.txt \
--out /out_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.${OUTNAME}

cat ${outdir}/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.${OUTNAME}.fam | cut -d"_" -f1,2 | awk ' { print $1,"unk",$2 }' > ${outdir}/indfile.txt

numk=${NUMK}

cd $outdir

for K in $(seq 1 $numk); do
    CMD="singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir ${gwas_tools_image} bash -c '/admixture/admixture --cv -j8 /out_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.${OUTNAME}.bed $K | tee /out_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.${OUTNAME}.${K}.out'" # Normally you should give --cv as first option to admixture
    echo sbatch --time=04:00:00 -c 8 -J ADMIXTURE_${K} -A data-machine -o ${root}"/output_data/slurm_logs/13_Structure/popgen_ADMIXTURE_${K}.log" --export=root=${root} --mem 12000 --wrap="$CMD"
    sbatch --time=04:00:00 -c 8 -J ADMIXTURE_${K} -A data-machine -o ${root}"/output_data/slurm_logs/13_Structure/popgen_ADMIXTURE_${K}.log" --export=root=${root} --mem 12000 --wrap="$CMD"
done
