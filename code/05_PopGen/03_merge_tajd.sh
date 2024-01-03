## 03_merge_tajd.sh
## 12-03-23
## JRG
## Merges chromsome split Tajima's D Files into Population TajD Files

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

POPS=($(cat ${root}/input_data/05_PopGen/popsfile.txt | cut -f2 | uniq))

outdir=${root}/output_data/05_PopGen/

for p in "${POPS[@]}"
do
   #tajd
   echo "Working on population $p..."
   cat ${outdir}/tajd_chunks/${p}_TAJD* > ${outdir}/tajd.${p}.combined.txt
   header=$(head -1 ${outdir}/tajd.${p}.combined.txt)
   sed -i "s/$header//g" ${outdir}/tajd.${p}.combined.txt
   echo -e $header | cat - ${outdir}/tajd.${p}.combined.txt > ${outdir}/${p}.tajd.stats.txt
   rm ${outdir}/tajd.${p}.combined.txt
   echo "...Done!"
done