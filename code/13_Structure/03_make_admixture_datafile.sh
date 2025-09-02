#!/bin/bash

## 07_make_admixture_datafile.sh
## Jan 2024
## JRG
# This script combines the ADMIXTURE datafiles for the different K clusters

root="$(git rev-parse --show-toplevel)"

source $root/pkings_gwas.env

MYNAME=$(whoami)

OUTNAME=$1

indir=${root}/output_data/13_Structure/$OUTNAME/
outdir=${root}/output_data/13_Structure/$OUTNAME/

mkdir -p ${outdir}

OUTROOT=${snps_only_vcf%".vcf.gz"}

rm -f ${outdir}/cv_errors.txt
grep -h CV ${outdir}/${OUTROOT}.renamed.maf5.miss20.dp5.r20kb.for_admixture.$OUTNAME.*.out | cut -d"=" -f2 | sed --expression='s/):/\t/g' > ${outdir}/cv_errors.txt


rm -f ${outdir}/admixture_long_results.txt
for i in $(cut -f1 ${outdir}/cv_errors.txt);
  do
   for k in $(seq 1 $i);
	do
      cat  ${outdir}/${OUTROOT}.renamed.maf5.miss20.dp5.r20kb.for_admixture.$OUTNAME.${i}.Q | cut -d" " -f$k | awk -v i="$i" -v k="$k" '{print i,k,$0}'| paste ${outdir}/indfile.txt - >> ${outdir}/admixture_long_results.txt
   done
done

echo -e "Sample\tSex\tPop\tMAXK\tGROUP\tQ" | cat - ${outdir}/admixture_long_results.txt > ${outdir}/admixture_long_results.header.txt
sed -i 's/ /\t/g' ${outdir}/admixture_long_results.header.txt