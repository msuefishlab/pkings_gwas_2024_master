#!/bin/bash --login

OUTNAME=$1

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

outdir=${root}/output_data/07_Peak_Analysis/${OUTNAME}

mkdir -p ${outdir}

# Rscript ${scriptdir}/render_report.R INPUTDATA OUTPUTDATA REPORT.Rmd
singularity exec ${gwas_tools_image} Rscript ${root}/code/00_common/render_rmd_report.R ${OUTNAME} ${outdir}/${OUTNAME}_Peak_Analysis ${root}/code/08_Peak_Analysis/peak_analysis.Rmd
