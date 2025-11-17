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
- Configuration file `tfbs.env` with paths to reference genome and GFF annotation

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
Each run creates a directory `output_data/<FRIENDLY_NAME>/` containing:
- `*_region.bed`: Genomic region analyzed
- `*_BP.fa` / `*_TP.fa`: Extracted sequences
- `*_fimo_compare/`: FIMO results (TSV, HTML, etc.)
- `*_motifs_BP_TP.bed`: All motifs with haplotype labels
- `*_motif_SNP_BP_TP_overlap.bed`: Motifs overlapping SNPs
- `*_motif_snp_summary.tsv`: Summary table
- `*_motif_snp_{BP_only,TP_only,both}.tsv`: Categorized results

### Python Analysis Script

The `analyze_motif_snp_overlap.py` script:
- Takes bedtools intersect output (motifs × SNPs)
- Aggregates best FIMO score per SNP-motif-haplotype combination
- Pivots to compare BP vs TP scores
- Outputs summary statistics and categorized subsets

## Environment Configuration

Edit `tfbs.env` to set:
- `SOURCE_GFF`: Path to genome annotation (GFF format)
- `REF`: Path to reference genome FASTA file

## Input Data Requirements

Required files in `input_data/`:
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
