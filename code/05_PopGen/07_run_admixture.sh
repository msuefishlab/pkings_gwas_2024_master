#!/bin/bash

## 03_run_admixture.sh
## Jan 2024
## JRG
# This script runs ADMIXTURE on the GDS file genereated for the phylogeny

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

MYNAME=$(whoami)

indir=${root}/output_data/05_PopGen/
keepdir=${root}/output_data/05_PopGen/
outdir=${root}/output_data/05_PopGen/

mkdir -p ${outdir}

OUTROOT=${snps_only_vcf%".vcf.gz"}

## keep only P. kingsleyae
singularity exec --bind $root:/project_root --bind $keepdir:/keep_dir --bind $outdir:/out_dir --bind $indir:/in_dir ${gwas_tools_image} /plink/plink \
--bfile /out_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture  \
--allow-extra-chr \
--allow-no-sex \
--make-bed \
--keep /keep_dir/pkings_keep.keep.txt \
--out /out_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.kings_only

cat ${outdir}/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.kings_only.fam | cut -d"_" -f1,2 | awk ' { print $1,"unk",$2 }' > ${outdir}/indfile.txt

numk=$(cat ${outdir}/indfile.txt | cut -d" " -f3 | uniq | wc -l)

cd $outdir

for K in $(seq 1 $numk); do
    CMD="singularity exec --bind $root:/project_root --bind $outdir:/out_dir --bind $indir:/in_dir ${gwas_tools_image} bash -c '/admixture/admixture --cv -j8 /out_dir/$OUTNAME/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.kings_only.bed $K | tee /out_dir/$OUTROOT.renamed.maf5.miss20.dp5.r20kb.for_admixture.kings_only.${K}.out'" # Normally you should give --cv as first option to admixture
    echo sbatch --time=04:00:00 -c 8 -J ADMIXTURE_${K} -o ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/popgen_ADMIXTURE_${K}.log" --export=root=${root} --mem 12000 --wrap="$CMD"
    sbatch --time=04:00:00 -c 8 -J ADMIXTURE_${K} -o ${root}"/output_data/slurm_logs/04_Phylogeny_And_Admixture/popgen_ADMIXTURE_${K}.log" --export=root=${root} --mem 12000 --wrap="$CMD"
done
