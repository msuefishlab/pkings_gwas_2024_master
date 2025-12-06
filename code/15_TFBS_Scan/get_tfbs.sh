#!/bin/bash
## get tfbs.sh
## 11-14-2025
## JRG
## Gets TFBS from BP and TP Individuals

set -eo pipefail  # Exit on error and pipe failures (but allow unset variables for now)

root="$(git rev-parse --show-toplevel)"
source "${root}/pkings_gwas.env"

set -u  # Now enable strict variable checking

# Activate meme conda environment (try conda first, then mamba)
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate meme
elif command -v mamba &> /dev/null; then
    eval "$(mamba shell.bash hook)"
    mamba activate meme
else
    echo "Error: Neither conda nor mamba found. Please install conda/mamba." >&2
    exit 1
fi

GENE=$1
FRIENDLY_NAME=$2
SAMPLE=$3

# Validate inputs
if [[ -z "$GENE" || -z "$FRIENDLY_NAME" || -z "$SAMPLE" ]]; then
    echo "Error: Missing required arguments" >&2
    echo "Usage: $0 <GENE_OR_REGION> <FRIENDLY_NAME> <SAMPLE_ID>" >&2
    exit 1
fi

mkdir -p output_data/15_TFBS_Scan/${FRIENDLY_NAME}

# Check if GENE is a genomic region (chr:start-end) or a gene name
if [[ $GENE =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
    # Input is a genomic region
    CHR=${BASH_REMATCH[1]}
    REGION_START=${BASH_REMATCH[2]}  # 1-based
    REGION_END=${BASH_REMATCH[3]}    # 1-based

    # Validate coordinates
    if (( REGION_START >= REGION_END )); then
        echo "Error: Invalid coordinates. Start ($REGION_START) must be less than end ($REGION_END)" >&2
        echo "Note: Coordinates should be in format chr:start-end (e.g., chr8:1000000-2000000)" >&2
        exit 1
    fi

    # Create BED file (0-based start, 1-based end)
    printf "%s\t%d\t%d\t%s\n" "${CHR}" "$((REGION_START-1))" "${REGION_END}" "${FRIENDLY_NAME}" \
        > output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_region.bed

    # Set variables for downstream use
    START=$((REGION_START-1))  # 0-based for BED
    END=${REGION_END}
else
    # Input is a gene name - extract from GFF
    grep -w ${GENE} ${TFBS_GFF} > output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}.gff

    # Get min start / max end
    awk -v fname="${FRIENDLY_NAME}" '
    NR==1 { chr=$1; start=$4; end=$5; next }
    { if ($4 < start) start=$4; if ($5 > end) end=$5 }
    END { printf "%s\t%d\t%d\t%s\n", chr, start-1, end, fname }
    ' output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}.gff > output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_region.bed

    read CHR START END FRIENDLY_NAME < output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_region.bed
fi

bedtools getfasta \
-fi ${TFBS_REF} \
-bed output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_region.bed \
-fo output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.fa

bcftools view -r ${CHR}:${START}-${END} -Oz -o output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.vcf.gz ${TFBS_VCF}
bcftools index output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.vcf.gz

bcftools view -r ${CHR}:${START}-${END} -s ${SAMPLE} -Oz -o output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.vcf.gz ${TFBS_VCF}
bcftools index output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.vcf.gz

samtools faidx ${TFBS_REF} ${CHR}:${START}-${END} | bcftools consensus -s ${SAMPLE} output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.vcf.gz > output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.fa

# Make a BED of SNPs with p-values

bcftools query -f '%CHROM\t%POS\n' output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.vcf.gz \
| awk 'NR==FNR {log_p[$1":"$3]=$13; rs_id[$1":"$3]=$2; next}
{key=$1":"$2; print $1"\t"($2-1)"\t"$2"\t"(rs_id[key] ? rs_id[key] : ".")"\t"(log_p[key] ? log_p[key] : "NA")}' \
 <(tail -n +2 ${TFBS_SNP_METADATA}) - \
> output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_snps.bed

BP_REGION=${CHR}:${START}-${END}
TP_REGION=${CHR}:${START}-${END}

sed "1s/^>.*/>BP|$BP_REGION/" output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.fa > output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.fa.tmp && mv output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.fa.tmp output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.fa
sed "1s/^>.*/>TP|$TP_REGION/" output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.fa > output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.fa.tmp && mv output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.fa.tmp output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.fa

cat output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_BP.fa output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_TP.fa > output_data/15_TFBS_Scan/${FRIENDLY_NAME}/bp_tp_${FRIENDLY_NAME}_cat.fa

fimo --oc output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_fimo_compare \
 --thresh 1e-4 \
 ${JASPAR_MOTIFS} \
 output_data/15_TFBS_Scan/${FRIENDLY_NAME}/bp_tp_${FRIENDLY_NAME}_cat.fa

awk '
BEGIN { OFS="\t"; header_seen=0 }
/^#/ { next } # skip comment lines
/^[[:space:]]*$/ { next } # skip blank lines
!header_seen { header_seen=1; next } # skip first non-comment, non-blank line (header)
{
motif = $2
seq = $3 # e.g. BP|chr17:11672113-11736094
start = $4 # 1-based within sequence
end = $5 # 1-based within sequence (inclusive)
score = $7
pval = $8

# split seq into hap and coordinates

n = split(seq, a, "|")
hap = a[1]
coord = a[2]

# split coordinates into chr and region span

split(coord, b, ":")
chr = b[1]
split(b[2], c, "-")
region_start = c[1] + 0 # 1-based genomic start of the sequence

# convert to genomic coordinates (1-based)

g_start_1 = region_start + start - 1
g_end_1 = region_start + end - 1

# BED: 0-based start, half-open end

bed_start = g_start_1 - 1
bed_end = g_end_1

# chr, start, end, motif, score, hap, pval

printf "%s\t%d\t%d\t%s\t%f\t%s\t%s\n", chr, bed_start, bed_end, motif, score, hap, pval
}' output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_fimo_compare/fimo.tsv \
> output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_motifs_BP_TP.bed



bedtools intersect -wa -wb \
 -a output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_motifs_BP_TP.bed \
 -b output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_snps.bed \
> output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_motif_SNP_BP_TP_overlap.bed

python3 code/15_TFBS_Scan/analyze_motif_snp_overlap.py \
 --input output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}_motif_SNP_BP_TP_overlap.bed \
 --output_prefix output_data/15_TFBS_Scan/${FRIENDLY_NAME}/${FRIENDLY_NAME}
