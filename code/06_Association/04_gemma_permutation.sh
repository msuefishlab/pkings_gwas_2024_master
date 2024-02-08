#!/bin/bash

## 03_gemma_permution.sh
## Feb 2023
## JRG
# This script submits a previously run association analysis (arg1) with phenotypes in arg2 to run a permuted GEMMA lmm asscoiation analysis
# e.g. bash 04_gemma_permution.sh APA_BAM_ONLY_BP1_TP0_WOBBLE9

# note to get slurm working with gemma wrapper, pay close attention to the versions and add these lines to the slurm command in the gemma-wrapper code:
# #SBATCH --cpus-per-task=4
# #SBATCH --mem=120g
# #SBATCH -o gemma_permute_%j.out

OUTNAME=$1

cd $root

source $root/pkings_gwas.env

MYNAME=$(whoami)

indir=${root}/output_data/06_Association/
outdir=${root}/output_data/06_Association/
tmpoutdir=$SCRATCH/pkings_gwas_2024_scratch/output_data/06_Assoc/

mkdir -p ${outdir}
mkdir -p ${tmpoutdir}

mkdir -p ${outdir}/$OUTNAME
mkdir -p ${tmpoutdir}/$OUTNAME

mkdir -p ${outdir}/${OUTNAME}_permution/

singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} bash -c "/gemma-wrapper/gemma-wrapper/bin/gemma-wrapper \
--no-parallel --cache-dir /out_dir/$OUTNAME/${OUTNAME}_permution \
--json -- \
-g /tmp_dir/$OUTNAME/${OUTNAME}.prune.for_gemma.geno \
-p /out_dir/$OUTNAME/$OUTNAME.pheno \
-gk 1 \
-debug > /out_dir/${OUTNAME}_permution/${OUTNAME}_K.json"


cd ${outdir}/${OUTNAME}_permution/

singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $tmpoutdir:/tmp_dir ${gwas_tools_image} bash -c "/gemma-wrapper/gemma-wrapper/bin/gemma-wrapper \
--input /out_dir/${OUTNAME}_permution/${OUTNAME}_K.json \
--permutate 100 --permute-phenotype /out_dir/$OUTNAME/$OUTNAME.pheno \
--cache-dir /out_dir/$OUTNAME/${OUTNAME}_permution \
--slurm -- \
-g /out_dir/$OUTNAME/$OUTNAME.geno \
-a /out_dir/$OUTNAME/${OUTNAME}_merge_map.txt \
-lmm 2 \
-debug > ${OUTNAME}_GWA_100_lmm2.txt"