#!/bin/bash --login

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

MYNAME=$(whoami)

outdir=${root}/output_data/04_Phylogenomics/

mkdir -p ${outdir}

OUTROOT=${snps_only_vcf%".vcf.gz"}

Rscript ${root}/code/00_common/render_rmd_report.R NULL ${outdir}/${OUTROOT}_mormyrid_phylogeny ${root}/code/04_Phylogenomics/phylogeny_and_character_reconstruction_mormyrids.Rmd

Rscript ${root}/code/00_common/render_rmd_report.R NULL ${outdir}/${OUTROOT}_paramormyrops_phylogeny ${root}/code/04_Phylogenomics/phylogeny_and_character_reconstruction_paramormyrops.Rmd