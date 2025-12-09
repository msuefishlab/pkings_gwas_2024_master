#!/bin/bash
set -euo pipefail

## polarize_vcf.sh
## August 2025
## JRG

module load BCFtools   # or conda env with bcftools≥1.13

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

outdir="${root}/output_data/12_Sweep_Detection"
mkdir -p "${outdir}"

vcf="${root}/output_data/07_Phasing/${phased_vcf}"
out_vcf="${outdir}/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.phased.renamed.polarized.vcf.gz"
outgroup="${outdir}/outgroup.txt"

echo "[1/5] Extracting outgroup sample IDs..."
bcftools query -l "${vcf}" | grep '^PSZA_' > "${outgroup}"
if [[ ! -s "${outgroup}" ]]; then
  echo "ERROR: outgroup list is empty (no samples starting with PSZA_ found)"; exit 1
fi

echo "[2/5] Counting alleles in the outgroup..."
bcftools view -S "${outgroup}" -Ou "${vcf}" \
  | bcftools view -g hom -Ou \
  | bcftools +fill-tags -Ou -- -t AN,AC \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\n' \
  > "${outdir}/outgroup_counts.tsv"

# If you prefer a looser rule (allow hom. but some missing), drop "-g hom" and enforce thresholds in awk below.

echo "[3/5] Deciding AA from outgroup counts..."
awk 'BEGIN{OFS="\t"}
{
  chrom=$1; pos=$2; ref=$3; alt=$4; an=$5; ac=$6;
  aa=".";
  # REQUIREMENTS (tune as desired):
  # - at least 8 outgroup chromosomes called
  # - and either AC==0 (fixed REF) or AC==AN (fixed ALT)
  if (an >= 8) {
    if (ac == 0)       aa=ref;
    else if (ac == an) aa=alt;
    else               aa=".";
  }
  print chrom, pos, aa
}' "${outdir}/outgroup_counts.tsv" > "${outdir}/AA.tsv"

# Compress & index the annotation table for bcftools annotate
bgzip -f "${outdir}/AA.tsv"
tabix -f -s 1 -b 2 -e 2 "${outdir}/AA.tsv.gz"

echo "[4/5] Writing INFO/AA onto the full VCF..."
bcftools annotate \
  -a "${outdir}/AA.tsv.gz" \
  -c CHROM,POS,INFO/AA \
  -h <(echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele inferred from PSZA_ outgroup (fixed across outgroup)">' ) \
  -Oz -o "${out_vcf}" "${vcf}"

echo "[5/5] Indexing result..."
bcftools index -f "${out_vcf}"

echo "Done. Polarized VCF: ${out_vcf}"
