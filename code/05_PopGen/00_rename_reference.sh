#!/bin/bash

## 00_rename_reference.sh
## Jan 2024
## JRG
# This script renames the reference fasta based on the chr_name_map for easier analysis

root="$(git rev-parse --show-toplevel)"
source ${root}/"pkings_gwas.env"

cat ${root}/input_data/00_Reference_Genome/${chr_name_map} | cut -f1,2 > ${root}/input_data/00_Reference_Genome/${chr_name_map}.only.txt
singularity exec --bind $root:/project_root ${gwas_tools_image}  /seqkit/seqkit replace -p '^(\S+)(.+?)$' -r '{kv}$2' -k /project_root/input_data/00_Reference_Genome/${chr_name_map}.only.txt /project_root/input_data/00_Reference_Genome/${reference} > /project_root/input_data/00_Reference_Genome/${reference}.renamed.fa