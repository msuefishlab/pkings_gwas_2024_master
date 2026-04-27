#!/bin/bash
## get_tfbs_background.sh
## Background window version of get_tfbs.sh for null distribution analysis.
##
## Differences from get_tfbs.sh:
##   1. Uses $phased_VCF (phased full-genome SNPs) instead of $TFBS_VCF
##   2. SNP BED is built from bcftools query with log_p=NA (no GWAS metadata)
##   3. Output goes to output_data/15_TFBS_Scan/background/<FRIENDLY_NAME>/
##
## Everything else is identical to get_tfbs.sh.
##
## Usage:
##   bash code/15_TFBS_Scan/get_tfbs_background.sh <REGION> <WINDOW_ID> <SAMPLE_ID>
##   e.g.: bash ... chr17:11000000-11100000 bg_sptbn4b_0001 APA_1234

set -eo pipefail

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

set -u

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate meme
elif command -v mamba &> /dev/null; then
    eval "$(mamba shell.bash hook)"
    mamba activate meme
else
    echo "Error: Neither conda nor mamba found." >&2
    exit 1
fi

GENE=$1
FRIENDLY_NAME=$2
SAMPLE=$3

if [[ -z "$GENE" || -z "$FRIENDLY_NAME" || -z "$SAMPLE" ]]; then
    echo "Usage: $0 <REGION> <WINDOW_ID> <SAMPLE_ID>" >&2
    exit 1
fi

OUTDIR="${root}/output_data/15_TFBS_Scan/background/${FRIENDLY_NAME}"
mkdir -p "${OUTDIR}"

# ── Region parsing (genomic coords only — background windows are always regions)
if [[ $GENE =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
    CHR=${BASH_REMATCH[1]}
    REGION_START=${BASH_REMATCH[2]}
    REGION_END=${BASH_REMATCH[3]}
else
    echo "Error: background windows must be genomic regions (chr:start-end), not gene names." >&2
    exit 1
fi

if (( REGION_START >= REGION_END )); then
    echo "Error: Invalid coordinates. Start ($REGION_START) must be less than end ($REGION_END)" >&2
    exit 1
fi

printf "%s\t%d\t%d\t%s\n" "${CHR}" "$((REGION_START-1))" "${REGION_END}" "${FRIENDLY_NAME}" \
    > "${OUTDIR}/${FRIENDLY_NAME}_region.bed"

START=$((REGION_START-1))
END=${REGION_END}

# ── Sequence extraction (identical to get_tfbs.sh) ────────────────────────────
bedtools getfasta \
    -fi "${TFBS_REF}" \
    -bed "${OUTDIR}/${FRIENDLY_NAME}_region.bed" \
    -fo "${OUTDIR}/${FRIENDLY_NAME}_BP.fa"

# ── Variant extraction — use full-genome SNP VCF, not TFBS_VCF ───────────────
bcftools view -r "${CHR}:${START}-${END}" -Oz \
    -o "${OUTDIR}/${FRIENDLY_NAME}_BP.vcf.gz" "${phased_VCF}"
bcftools index "${OUTDIR}/${FRIENDLY_NAME}_BP.vcf.gz"

bcftools view -r "${CHR}:${START}-${END}" -s "${SAMPLE}" -Oz \
    -o "${OUTDIR}/${FRIENDLY_NAME}_TP.vcf.gz" "${phased_VCF}"
bcftools index "${OUTDIR}/${FRIENDLY_NAME}_TP.vcf.gz"

samtools faidx "${TFBS_REF}" "${CHR}:${START}-${END}" \
    | bcftools consensus -s "${SAMPLE}" "${OUTDIR}/${FRIENDLY_NAME}_TP.vcf.gz" \
    > "${OUTDIR}/${FRIENDLY_NAME}_TP.fa"

# ── SNP BED — positions only, log_p=NA (no GWAS metadata for background SNPs)
bcftools query -r "${CHR}:${START}-${END}" \
    -f '%CHROM\t%POS\n' \
    "${OUTDIR}/${FRIENDLY_NAME}_BP.vcf.gz" \
| awk '{print $1"\t"($2-1)"\t"$2"\t"$1"_"$2"\tNA"}' \
> "${OUTDIR}/${FRIENDLY_NAME}_snps.bed"

# ── Label FASTA headers (identical to get_tfbs.sh) ───────────────────────────
BP_REGION="${CHR}:${START}-${END}"
TP_REGION="${CHR}:${START}-${END}"

sed "1s/^>.*/>BP|$BP_REGION/" "${OUTDIR}/${FRIENDLY_NAME}_BP.fa" \
    > "${OUTDIR}/${FRIENDLY_NAME}_BP.fa.tmp" \
    && mv "${OUTDIR}/${FRIENDLY_NAME}_BP.fa.tmp" "${OUTDIR}/${FRIENDLY_NAME}_BP.fa"

sed "1s/^>.*/>TP|$TP_REGION/" "${OUTDIR}/${FRIENDLY_NAME}_TP.fa" \
    > "${OUTDIR}/${FRIENDLY_NAME}_TP.fa.tmp" \
    && mv "${OUTDIR}/${FRIENDLY_NAME}_TP.fa.tmp" "${OUTDIR}/${FRIENDLY_NAME}_TP.fa"

cat "${OUTDIR}/${FRIENDLY_NAME}_BP.fa" "${OUTDIR}/${FRIENDLY_NAME}_TP.fa" \
    > "${OUTDIR}/bp_tp_${FRIENDLY_NAME}_cat.fa"

# ── FIMO motif scanning (identical to get_tfbs.sh) ───────────────────────────
fimo --oc "${OUTDIR}/${FRIENDLY_NAME}_fimo_compare" \
     --thresh 1e-4 \
     "${JASPAR_MOTIFS}" \
     "${OUTDIR}/bp_tp_${FRIENDLY_NAME}_cat.fa"

# ── FIMO output → genomic BED (identical AWK to get_tfbs.sh) ─────────────────
awk '
BEGIN { OFS="\t"; header_seen=0 }
/^#/ { next }
/^[[:space:]]*$/ { next }
!header_seen { header_seen=1; next }
{
motif = $2
seq = $3
start = $4
end = $5
score = $7
pval = $8

n = split(seq, a, "|")
hap = a[1]
coord = a[2]

split(coord, b, ":")
chr = b[1]
split(b[2], c, "-")
region_start = c[1] + 0

g_start_1 = region_start + start - 1
g_end_1 = region_start + end - 1

bed_start = g_start_1 - 1
bed_end = g_end_1

printf "%s\t%d\t%d\t%s\t%f\t%s\t%s\n", chr, bed_start, bed_end, motif, score, hap, pval
}' "${OUTDIR}/${FRIENDLY_NAME}_fimo_compare/fimo.tsv" \
> "${OUTDIR}/${FRIENDLY_NAME}_motifs_BP_TP.bed"

# ── Intersect motifs with SNPs (identical to get_tfbs.sh) ────────────────────
bedtools intersect -wa -wb \
    -a "${OUTDIR}/${FRIENDLY_NAME}_motifs_BP_TP.bed" \
    -b "${OUTDIR}/${FRIENDLY_NAME}_snps.bed" \
> "${OUTDIR}/${FRIENDLY_NAME}_motif_SNP_BP_TP_overlap.bed"

# ── Python analysis (identical to get_tfbs.sh) ───────────────────────────────
python3 "${root}/code/15_TFBS_Scan/analyze_motif_snp_overlap.py" \
    --input "${OUTDIR}/${FRIENDLY_NAME}_motif_SNP_BP_TP_overlap.bed" \
    --output_prefix "${OUTDIR}/${FRIENDLY_NAME}"
