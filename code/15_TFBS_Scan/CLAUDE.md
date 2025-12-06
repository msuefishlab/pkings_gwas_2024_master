# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains a pipeline for analyzing transcription factor binding site (TFBS) differences between two haplotypes (BP and TP) at SNP locations. The pipeline identifies motifs that are gained or lost due to genetic variants.

## Running the Pipeline

### Main Command

```bash
bash code/get_tfbs.sh <GENE_OR_REGION> <FRIENDLY_NAME> <SAMPLE_ID>
```

**Input formats:**
- Gene name: `bash code/get_tfbs.sh SPTBN4 sptbn4 SAMPLE123`
- Genomic region (1-based coords): `bash code/get_tfbs.sh chr17:1000000-2000000 my_region SAMPLE123`

### Prerequisites

The pipeline requires:
- Conda environment named `meme` (activated automatically by the script)
- Bioinformatics tools: `bcftools`, `samtools`, `bedtools`, `fimo` (from MEME suite)
- Python 3 with `pandas`
- Environment configuration in `pkings_gwas.env` (sourced automatically by the script)

## Architecture

### Pipeline Flow (get_tfbs.sh)

1. **Region Extraction**: Parse gene name from GFF or use provided genomic coordinates
2. **Sequence Extraction**: Extract reference sequences for both BP and TP haplotypes
   - BP: Population-level VCF (all SNPs in peaks)
   - TP: Test individual's personal genome (consensus sequence)
3. **Motif Scanning**: Run FIMO (threshold: 1e-4) using JASPAR 2024 vertebrate motifs
4. **SNP Overlap**: Identify which motifs overlap with significant SNPs
5. **BP vs TP Analysis**: Python script categorizes motifs as BP-only, TP-only, or both

### Key Technical Details

**Coordinate Systems:**
- Input: 1-based genomic coordinates
- BED files: 0-based start, 1-based end (standard BED format)
- The script handles conversion automatically

**Haplotype Comparison:**
- BP (base population): Reference genome + population VCF
- TP (test population/individual): Sample-specific consensus sequence
- FASTA headers labeled as `>BP|chr:start-end` and `>TP|chr:start-end`

**Output Structure:**
Each run creates a directory `output_data/15_TFBS_Scan/<FRIENDLY_NAME>/` containing:
- `*_region.bed`: Genomic region analyzed
- `*_BP.fa` / `*_TP.fa`: Extracted sequences
- `*_fimo_compare/`: FIMO results (TSV, HTML, etc.)
- `*_motifs_BP_TP.bed`: All motifs with haplotype labels
- `*_motif_SNP_BP_TP_overlap.bed`: Motifs overlapping SNPs
- `*_motif_snp_summary.tsv`: Summary table
- `*_motif_snp_{BP_only,TP_only,both}.tsv`: Categorized results

### Python Analysis Script

The `analyze_motif_snp_overlap.py` script is the final step of the pipeline, processing the intersection of motifs and SNPs to categorize haplotype-specific binding site differences.

**Location:** `code/15_TFBS_Scan/analyze_motif_snp_overlap.py`

**Command Line:**
```bash
python code/15_TFBS_Scan/analyze_motif_snp_overlap.py \
  --input <bedtools_intersect_output.bed> \
  --output_prefix <output_directory/prefix>
```

**Required Arguments:**
- `--input` / `-i`: Path to 12-column BED file from `bedtools intersect -wa -wb`
- `--output_prefix` / `-o`: Output path prefix (e.g., `output_data/sptbn4b/sptbn4b`)

**Input Format (12 columns from bedtools intersect):**
- Columns 0-6 (from motifs_BP_TP.bed): chr, start, end, motif_id, fimo_score, haplotype, fimo_pval
- Columns 7-11 (from snps.bed): chr, start, end, snp_id, snp_log10p

**Processing Logic:**
1. **Aggregation**: Groups by (snp_id, motif_id, haplotype) and takes max FIMO score
   - Handles cases where a motif may overlap a SNP multiple times
2. **Pivoting**: Reshapes data to compare BP vs TP side-by-side
   - Creates separate columns: BP_score, TP_score, BP_present, TP_present
3. **Score Delta**: Calculates `score_delta_TP_minus_BP` for motifs present in both
   - Positive delta = stronger binding in TP
   - Negative delta = stronger binding in BP
4. **Sorting**: Orders by SNP significance (snp_log10p descending), then by snp_id and motif_id

**Output Files:**
- `<prefix>_motif_snp_summary.tsv`: All SNP-motif pairs with comparative scores
- `<prefix>_motif_snp_TP_only.tsv`: Motifs gained in TP (TP_present=True, BP_present=False)
- `<prefix>_motif_snp_TP_only.bed`: BED file of TP-specific motif locations
- `<prefix>_motif_snp_BP_only.tsv`: Motifs lost in TP (BP_present=True, TP_present=False)
- `<prefix>_motif_snp_BP_only.bed`: BED file of BP-specific motif locations
- `<prefix>_motif_snp_both.tsv`: Motifs in both (BP_present=True, TP_present=True)
- `<prefix>_motif_snp_both.bed`: BED file of shared motif locations (scored by SNP significance)

**BED File Details:**
- Standard 6-column BED format (chr, start, end, name, score, strand)
- Coordinates are motif locations (0-based BED format)
- Name field: `{motif_id}_{snp_id}_{haplotype}`
- Score field: FIMO score × 100 (for BP/TP-only files) or SNP log10p × 100 (for both file)
- Strand: always "." (FIMO scans both strands)
- Useful for IGV visualization, bedtools operations, or genomic overlap analyses

**Dependencies:**
- Python 3 with pandas
- No other special requirements

**Common Issues:**
- If input file is empty, script exits with error code 1
- Missing haplotype labels (BP/TP) in input will cause KeyError - check FASTA headers
- Non-numeric log10p values are coerced to NaN (won't cause failure)

## Environment Configuration

The pipeline uses variables defined in `pkings_gwas.env`:
- `TFBS_GFF`: Genome annotation (GFF format)
- `TFBS_REF`: Reference genome FASTA
- `TFBS_VCF`: VCF file with SNPs in peaks
- `TFBS_SNP_METADATA`: SNP metadata with p-values
- `JASPAR_MOTIFS`: JASPAR motif database (MEME format)

These are automatically sourced by `get_tfbs.sh` - no manual configuration needed.

## Input Data Requirements

Required files in `input_data/15_TFBS_Scan/`:
- `PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz`: Population VCF with SNPs
- `PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt`: SNP metadata with p-values (tab-delimited, header row)
- `JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt`: MEME-formatted motif database

## Modifying the Pipeline

**Bash Script (`get_tfbs.sh`):**
- AWK scripts handle GFF parsing and FIMO output processing
- Comment/blank line filtering is critical for FIMO TSV parsing
- Temporary files (`.tmp`) are used for in-place sed operations to avoid data loss

**FIMO Processing:**
- AWK script skips comment lines (`/^#/`), blank lines, and header
- Converts FIMO 1-based sequence coordinates to genomic BED coordinates
- Preserves haplotype labels (BP/TP) from FASTA headers
