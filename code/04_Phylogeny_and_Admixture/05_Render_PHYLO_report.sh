#!/bin/bash --login

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

MYNAME=$(whoami)

outdir=${root}/output_data/04_Phylogeny_And_Admixture/

mkdir -p ${outdir}

OUTROOT=${snps_only_vcf%".vcf.gz"}

Rscript ${root}/code/00_common/render_rmd_report.R NULL ${outdir}/${OUTROOT}_phylogeny ${root}/code/04_Phylogeny_and_Admixture/phylogeny_and_character_reconstruction_chronos_ultrametric.Rmd