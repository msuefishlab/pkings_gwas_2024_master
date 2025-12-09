#!/bin/bash
set -euo pipefail

## split_vcfs_by_lists.sh
## August 2025
## JRG

module load BCFtools   # or conda env with bcftools ≥1.13

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

indir="${root}/input_data/12_Sweep_Detection"
outdir="${root}/output_data/12_Sweep_Detection"

# INPUTS
vcf="${outdir}/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.phased.renamed.polarized.vcf.gz"  # polarized VCF with INFO/AA present
exclude_regex='^$'  # optional: e.g., '^PSZA$|^COB$' or leave empty '' / '^$' for no exclusion

mkdir -p "${outdir}"/{lists,vcfs}

# cache VCF sample names (sorted, unique)
vcf_samples_cache="${outdir}/lists/_ALL_VCF_SAMPLES.txt"
bcftools query -l "${vcf}" | sort -u > "${vcf_samples_cache}"

echo "Reading sample sets from: ${indir}"
shopt -s nullglob
list_files=( "${indir}"/*.txt )
shopt -u nullglob

if (( ${#list_files[@]} == 0 )); then
  echo "ERROR: No .txt files found in ${indir}" >&2
  exit 1
fi

echo "Found ${#list_files[@]} set file(s):"
printf ' - %s\n' "${list_files[@]}"

for f in "${list_files[@]}"; do
  setname="$(basename "${f}" .txt)"

  # Normalize the provided list:
  # - keep first field per line (defensive), strip CRs/blank lines
  # - unique, then apply optional exclusion
  # - intersect with actual VCF samples so we never pass unknown IDs to bcftools
  tmp_clean="$(mktemp)"
  awk '{print $1}' "${f}" | sed 's/\r$//' | grep -v '^[[:space:]]*$' | sort -u > "${tmp_clean}"

  if [[ -n "${exclude_regex}" ]]; then
    grep -Ev "${exclude_regex}" "${tmp_clean}" > "${tmp_clean}.keep" || true
    mv "${tmp_clean}.keep" "${tmp_clean}"
  fi

  listfile="${outdir}/lists/${setname}.samples.txt"
  # Intersection (names present in the VCF)
  comm -12 "${vcf_samples_cache}" "${tmp_clean}" > "${listfile}" || true
  rm -f "${tmp_clean}"

  ns=$(wc -l < "${listfile}" | tr -d '[:space:]')
  if (( ns == 0 )); then
    echo "WARN: ${setname} has 0 matching samples in VCF—skipping."
    rm -f "${listfile}"
    continue
  fi

  echo "Splitting ${setname} (n=${ns}) ..."
  outvcf="${outdir}/vcfs/${setname}.polarized.vcf.gz"

  # Subset and keep polarized INFO; ensure biallelic SNPs; fill AC/AN; drop sites with AN==0
  bcftools view -S "${listfile}" -Ou "${vcf}" \
  | bcftools view -m2 -M2 -v snps -Ou \
  | bcftools +fill-tags -Ou -- -t AC,AN \
  | bcftools view -i 'INFO/AN>0' -Oz -o "${outvcf}"

  bcftools index -f "${outvcf}"
done

echo "Done. Outputs in ${outdir}/vcfs/"
