# TFBS Scan Pipeline

A pipeline for identifying transcription factor binding site (TFBS) differences between haplotypes at SNP locations.

## Overview

This pipeline compares transcription factor binding motifs between two haplotypes (BP and TP) to identify motifs that are gained or lost due to genetic variants. It extracts genomic sequences, scans for TF motifs using FIMO, and analyzes differences at SNP locations.

## Prerequisites

- **Conda/Mamba** package manager
- **Bioinformatics tools**: `bcftools`, `samtools`, `bedtools`
- **Python 3** with `pandas`
- **Reference data**: genome FASTA and GFF annotation (configured in `tfbs.env`)

## Installation

### 1. Create the MEME Environment

```bash
mamba create -n meme -c bioconda -c conda-forge meme
```

This creates a conda environment named `meme` with the MEME suite (including FIMO) installed.

### 2. Get the JASPAR Motif Database

Download the JASPAR vertebrate motif database in MEME format (one-time setup):

1. Go to the [JASPAR website](https://jaspar.elixir.no/)
2. Download the vertebrate motif set in MEME format (e.g., `JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt`)
3. Place the downloaded file in the `input_data/` directory:
   ```bash
   cp JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt input_data/
   ```

## Setup

1. Configure paths in `tfbs.env`:
   ```bash
   export SOURCE_GFF=/path/to/genome.gff
   export REF=/path/to/reference.fasta
   ```

2. Ensure required input data files exist in `input_data/`:
   - `PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz` (VCF with SNPs)
   - `PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt` (SNP metadata with p-values)
   - `JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt` (motif database)

## Usage

```bash
bash code/get_tfbs.sh <GENE_OR_REGION> <FRIENDLY_NAME> <SAMPLE_ID>
```

**Parameters:**
- `GENE_OR_REGION`: Either a gene name or genomic region (chr:start-end, 1-based coordinates)
- `FRIENDLY_NAME`: Output directory name (created under `output_data/`)
- `SAMPLE_ID`: Sample identifier in the VCF file

### Example 1: Using a Gene Name

```bash
bash code/get_tfbs.sh sptbn4b sptbn4b APA_193
```

This will:
- Search for gene `sptbn4b` in the GFF file
- Extract the genomic region spanning all exons/features of that gene
- Analyze TFBS differences for sample `APA_193`
- Create output in `output_data/sptbn4b/`

### Example 2: Using a Genomic Region

```bash
bash code/get_tfbs.sh chr16:10728257-11082159 cntnap5 APA_193
```

This will:
- Use the specified genomic region directly (chromosome 16, positions 10728257-11082159)
- Analyze TFBS differences for sample `APA_193`
- Create output in `output_data/cntnap5/`

## Output Files

Each run creates a directory `output_data/<FRIENDLY_NAME>/` containing:

### Intermediate Files
- `*_region.bed` - Genomic region being analyzed
- `*_BP.fa` / `*_TP.fa` - Extracted sequences for BP and TP haplotypes
- `*_BP.vcf.gz` / `*_TP.vcf.gz` - Filtered VCF files for the region
- `*_snps.bed` - SNPs in the region with p-values
- `*_fimo_compare/` - Complete FIMO output (TSV, HTML reports)

### Final Results
- `*_motif_snp_summary.tsv` - Complete summary of all SNP-motif pairs
- `*_motif_snp_BP_only.tsv` - Motifs present only in BP haplotype
- `*_motif_snp_TP_only.tsv` - Motifs present only in TP haplotype
- `*_motif_snp_both.tsv` - Motifs present in both haplotypes (with score differences)

## Understanding the Results

The pipeline identifies three categories of motifs:

1. **BP-only**: Motifs disrupted in the test individual (TP) - potential loss of TF binding
2. **TP-only**: Motifs created in the test individual (TP) - potential gain of TF binding
3. **Both**: Motifs present in both, but potentially with different binding scores

The summary files are sorted by SNP significance (p-value), making it easy to prioritize variants with the strongest association signals.

## Pipeline Steps

1. Extract genomic region (from gene name or coordinates)
2. Get reference sequence for BP haplotype
3. Generate TP haplotype consensus sequence for the specified sample
4. Scan both sequences for TF motifs using FIMO (JASPAR 2024, threshold 1e-4)
5. Identify SNPs overlapping with predicted motifs
6. Compare BP vs TP motif occurrences and scores
7. Generate summary reports

## Notes

- Input coordinates should be 1-based (standard genome browser format)
- The pipeline uses the MEME suite's FIMO tool with a p-value threshold of 1e-4
- Motifs are scanned using the JASPAR 2024 vertebrate motif database
- BP = "base population" (reference + population variants)
- TP = "test population/individual" (sample-specific genome)
