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

module load GCC/12.2.0 OpenBLAS/0.3.21 GSL/2.7

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

mkdir -p ${outdir}/$OUTNAME/${OUTNAME}_permution/

export GEMMA_COMMAND="singularity exec --bind $PWD:$PWD ${gwas_tools_image} /gemma/GEMMA/bin/gemma" 

echo "making relmat..."
~/gemma-wrapper-0.99.1/bin/gemma-wrapper \
--no-parallel --cache-dir ${outdir}/$OUTNAME/${OUTNAME}_permution \
--json -- \
-a ${outdir}/$OUTNAME/${OUTNAME}_merge_map.txt \
-g ${tmpoutdir}/$OUTNAME/${OUTNAME}.prune.for_gemma.geno \
-p ${outdir}/$OUTNAME/$OUTNAME.pheno \
-gk 1 \
-debug > ${outdir}/$OUTNAME/${OUTNAME}_permution/${OUTNAME}_K.json


echo "starting permution..."
cd ${outdir}/$OUTNAME/${OUTNAME}_permution/

~/gemma-wrapper-0.99.1/bin/gemma-wrapper \
--input ${outdir}/$OUTNAME/${OUTNAME}_permution/${OUTNAME}_K.json \
--no-parallel --permutate 100 \
--permute-phenotypes ${outdir}/$OUTNAME/$OUTNAME.pheno \
--cache-dir ${outdir}/$OUTNAME/${OUTNAME}_permution \
--slurm \
-- \
-g ${outdir}/$OUTNAME/$OUTNAME.geno \
-a ${outdir}/$OUTNAME/${OUTNAME}_merge_map.txt \
-lmm 2 \
-debug > ${OUTNAME}_GWA_100_lmm2.txt
