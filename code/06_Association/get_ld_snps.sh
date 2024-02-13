#!/bin/bash

##usage: get_ld_snps.sh bfile keepfile ld-window ld-window-kb ld-window-r2 ld-snp-list outfile

module purge; module load GCC/7.3.0-2.30  OpenMPI/3.1.1 GEMMA/0.98.1 PLINK/1.9b_4.1-x86_64

plink \
--allow-extra-chr \
--allow-no-sex \
--bfile $1 \
--keep $2 \
--out $7 \
--r2 with-freqs \
--ld-window $3 \
--ld-window-kb $4 \
--ld-snp-list $6 \
--ld-window-r2 $5 \
--threads 2
