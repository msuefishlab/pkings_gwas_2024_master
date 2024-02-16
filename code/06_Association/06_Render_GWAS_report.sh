#!/bin/bash --login

OUTNAME=$1

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

outdir=${root}/output_data/06_Association/${OUTNAME}

mkdir -p ${outdir}

# Rscript ${scriptdir}/render_report.R INPUTDATA OUTPUTDATA REPORT.Rmd
singularity exec ${gwas_tools_image} Rscript ${root}/code/00_common/render_rmd_report.R ${OUTNAME} ${outdir}/${OUTNAME}_GWAS_Stage1 ${root}/code/06_Association/gwas_report.Rmd