# EHH Quick Start Guide

## Overview

This workflow performs **genome-wide Extended Haplotype Homozygosity (EHH) analysis** to determine how "exceptional" GWAS peak selection signals are relative to the rest of the genome. It calculates percentile rankings and empirical p-values for peak statistics.

**Populations:** BP1, TP1, BP2, TP2, BP3 (5 populations)
**Comparisons:** BP1_TP1, BP2_TP2, BP3_TP2 (3 cross-population comparisons)
**GWAS Peaks:** 6 peaks (chr6, chr8, chr13, chr16, chr17, chr24)

**Statistics computed:**
- **inES** - Integrated Extended Haplotype Statistic (raw EHH measure per population)
- **iHS** - Integrated Haplotype Score (standardized inES, mean=0, SD=1)
- **Rsb** - Relative Selection Between populations (cross-population comparison)

---

## Prerequisites

Before running genome-wide EHH analysis, you must have:

1. **Polarized VCFs** - Created by `polarize_vcf.sh`
2. **Per-population VCF splits** - Created by `split_vcfs.sh`
3. **GWAS peaks BED file** - From GWAS analysis (06_Association)

### One-Time Setup: Prepare VCFs

If you haven't already run the peak-specific EHH analysis, prepare the VCFs:

```bash
source pkings_gwas.env

# 1. Polarize VCF with ancestral allele (AA) field
bash code/12_Sweep_Detection/ehh/polarize_vcf.sh

# 2. Split into per-population VCFs
bash code/12_Sweep_Detection/ehh/split_vcfs.sh

# Optional: Split by chromosome for faster processing
bash code/12_Sweep_Detection/ehh/split_vcf_by_chrom.sh
```

**Output:**
- `output_data/12_Sweep_Detection/*.polarized.vcf.gz` (5 populations)
- `output_data/12_Sweep_Detection/by_chrom/*.polarized.chr*.vcf.gz` (optional)

---

## Running the Genome-Wide Analysis

### 1. Submit inES/iHS Scan Jobs (~2-3 hours wall time)

Compute genome-wide inES and iHS for all populations:

```bash
bash code/12_Sweep_Detection/ehh/01_submit_ines_scans.sh

# For full genome (optional, default is GWAS peak chromosomes only):
# bash code/12_Sweep_Detection/ehh/01_submit_ines_scans.sh --all-chromosomes

# Monitor jobs:
squeue -u $USER
```

**Jobs submitted:**
- GWAS peak chromosomes only: 30 jobs (5 populations × 6 chromosomes)
- Full genome: 125 jobs (5 populations × 25 chromosomes)

**Resources per job:**
- Time: 2 hours
- Memory: 12 GB
- CPUs: 1

**Output:** `output_data/12_Sweep_Detection/ehh/ines_scans/{POP}.{CHR}.ines.txt.gz`

**What it does:**
- Loads polarized VCF for one population, one chromosome
- Filters to MAF ≥ 0.05
- Computes inES using `rehh::scan_hh()`
- Standardizes to iHS (mean=0, SD=1)
- Saves compressed TSV with columns: CHR, POSITION, inES, iHS

---

### 2. Submit Rsb Scan Jobs (~1-2 hours wall time)

Compute genome-wide Rsb between population pairs:

```bash
# Wait for step 1 to complete, then:
bash code/12_Sweep_Detection/ehh/02_submit_rsb_scans.sh

# For full genome (optional):
# bash code/12_Sweep_Detection/ehh/02_submit_rsb_scans.sh --all-chromosomes

# Monitor jobs:
squeue -u $USER
```

**Jobs submitted:**
- GWAS peak chromosomes only: 18 jobs (3 comparisons × 6 chromosomes)
- Full genome: 75 jobs (3 comparisons × 25 chromosomes)

**Resources per job:**
- Time: 1 hour
- Memory: 8 GB
- CPUs: 1

**Output:** `output_data/12_Sweep_Detection/ehh/rsb_scans/{POP1}_{POP2}.{CHR}.rsb.txt.gz`

**What it does:**
- Loads pre-computed inES files for two populations
- Computes Rsb using `rehh::ines2rsb()` (following Tang et al. 2007)
- Saves compressed TSV with columns: CHR, POSITION, RSB_{POP1}_{POP2}

---

### 3. Merge Per-Chromosome Results (~30 minutes)

Combine per-chromosome files into genome-wide distributions:

```bash
# After all jobs complete:
bash code/12_Sweep_Detection/ehh/03_merge_ehh_results.sh
```

**Output:**
- `output_data/12_Sweep_Detection/ehh/merged/{POP}.genome_wide_ines.txt.gz` (5 files)
- `output_data/12_Sweep_Detection/ehh/merged/{COMP}.genome_wide_rsb.txt.gz` (3 files)
- `output_data/12_Sweep_Detection/ehh/merged/merge.log`

**What it does:**
- Concatenates all chromosomes per population/comparison
- Preserves headers from first file
- Reports total lines and file sizes
- Validates all files created successfully

---

### 4. Extract Peak Values & Calculate Percentiles (~15 minutes)

This single script performs three tasks:
1. Calculate genome-wide summary statistics (mean, SD, percentiles)
2. Extract EHH values within GWAS peaks
3. Calculate percentile rankings and empirical p-values

```bash
Rscript code/12_Sweep_Detection/ehh/04_extract_peak_ehh.R
```

**Outputs:**

*Genome-wide summaries:*
- `genome_wide_ines_summary.txt` - Percentiles (1st, 5th, 10th, 50th, 90th, 95th, 99th) for inES/iHS
- `genome_wide_rsb_summary.txt` - Percentiles for Rsb

*Peak extractions:*
- `peak_ines_values.txt` - All SNPs within 6 GWAS peaks (inES/iHS)
- `peak_ines_summary.txt` - Summary per peak (mean, median, max, SD, position of max)
- `peak_rsb_values.txt` - All SNPs within peaks (Rsb)
- `peak_rsb_summary.txt` - Summary per peak

*Percentile rankings:*
- `peak_ines_percentiles.txt` - Percentile rank of each peak's max value
- `peak_rsb_percentiles.txt` - Percentile rank for Rsb

*Empirical p-values:*
- `peak_ines_empirical_pvalues.txt` - Statistical significance (upper tail, two-tailed)
- `peak_rsb_empirical_pvalues.txt` - Statistical significance

**What it does:**
- For each peak, calculates: percentile = `mean(genome_wide <= peak_max) × 100`
- For each peak, calculates: empirical p-value = `mean(genome_wide >= peak_max)`
- Ranks peaks by genome-wide rank (1 = highest value)

**Example interpretation:**
> "Peak 2 in BP1 has max |iHS| = 2.345, which is in the 99.87th percentile genome-wide (only 0.13% of SNPs have higher |iHS|). Empirical p-value = 0.0013 (genome-wide significant)."

---

### 5. Generate HTML Report (~20 minutes)

Create comprehensive HTML report with visualizations:

```bash
bash code/12_Sweep_Detection/ehh/05_render_ehh_report.sh
```

**Output:** `output_data/12_Sweep_Detection/ehh/ehh_analysis.html`

**What it includes:**
1. **Genome-wide summary statistics** - Mean, SD, percentiles for all populations/comparisons
2. **Percentile heatmaps** - Visual summary of peak rankings (rows=peaks, columns=populations/comparisons)
3. **Empirical p-values** - Statistical significance table
4. **Distribution comparison plots** - Genome-wide density (gray) with peak values overlaid (red rug)
5. **Manhattan plots** - Genome-wide landscape with GWAS peaks highlighted (blue), percentile thresholds (dashed lines)
6. **Peak summary tables** - Top-ranked peaks sorted by percentile

**Generated figures (saved as PDFs):**
- `ihs_percentile_heatmap.pdf` - Heatmap of iHS percentiles
- `rsb_percentile_heatmap.pdf` - Heatmap of Rsb percentiles
- `ihs_distribution_comparison.pdf` - Genome-wide vs peak distributions (iHS)
- `rsb_distribution_comparison.pdf` - Genome-wide vs peak distributions (Rsb)
- `ihs_manhattan.pdf` - Genome-wide Manhattan plot (iHS)
- `rsb_manhattan.pdf` - Genome-wide Manhattan plot (Rsb)

---

### 6. View Results

```bash
# Copy to your local machine or open in browser on HPCC
open output_data/12_Sweep_Detection/ehh/ehh_analysis.html
```

---

## Interpreting Results

### Percentile Rankings

- **>90%**: Peak is in top 10% genome-wide (moderate signal)
- **>95%**: Peak is in top 5% genome-wide (strong signal)
- **>99%**: Peak is in top 1% genome-wide (exceptional signal)

### Empirical P-Values

- **p < 0.05**: Significant at 5% level
- **p < 0.01**: Genome-wide significant at 1% level
- **p < 0.001**: Highly significant

### Statistics Interpretation

- **High iHS**: Recent positive selection **within** a population (hard sweep)
- **High Rsb**: Differential selection **between** populations (population-specific adaptation)
- **High inES**: Strong haplotype extension (raw measure, unnormalized)

### Biological Context

GWAS peaks with:
- High iHS + High Rsb → Strong population-specific sweep at EOD locus
- High iHS + Low Rsb → Shared selection across populations
- Low iHS + High Rsb → Old sweep with population-specific allele frequency differences

---

## Troubleshooting

**No inES files after step 1:**
```bash
# Check SLURM logs:
ls output_data/slurm_logs/12_Sweep_Detection/ehh/
tail output_data/slurm_logs/12_Sweep_Detection/ehh/ines_BP1_chr6.*.log

# Verify polarized VCFs exist:
ls output_data/12_Sweep_Detection/*.polarized.vcf.gz
```

**Rsb jobs skipped:**
```bash
# Check that inES files exist:
ls output_data/12_Sweep_Detection/ehh/ines_scans/BP1.chr*.ines.txt.gz

# If missing, re-run step 1
```

**Merge script reports missing files:**
```bash
# Check merge log:
cat output_data/12_Sweep_Detection/ehh/merged/merge.log

# Verify all chromosomes completed:
ls output_data/12_Sweep_Detection/ehh/ines_scans/ | grep BP1 | wc -l
# Should be 6 (GWAS peaks) or 25 (full genome)
```

**R Markdown rendering fails:**
```bash
# Check that all prerequisite files exist:
ls output_data/12_Sweep_Detection/ehh/merged/*.genome_wide_*.txt.gz
ls output_data/12_Sweep_Detection/ehh/peak_*.txt
ls output_data/12_Sweep_Detection/ehh/genome_wide_*_summary.txt

# If missing, re-run step 4
```

**Out of memory errors:**
```bash
# Step 1 (inES) needs 12 GB - if failing, request more:
# Edit ehh/01_run_ines_scan.sb: #SBATCH --mem=16GB

# Step 2 (Rsb) needs 8 GB - increase if needed:
# Edit ehh/02_run_rsb_scan.sb: #SBATCH --mem=12GB
```

---

## Total Time Estimate

- VCF preparation (one-time, if not done): 2-3 hours
- inES/iHS scanning (GWAS peaks): 2-3 hours wall time (parallelized)
- Rsb scanning (GWAS peaks): 1-2 hours wall time
- Merge + analysis + report: 1 hour

**Total:** ~4-6 hours (GWAS peak chromosomes)
**Total (full genome):** ~6-10 hours (all 25 chromosomes)

---

## File Sizes

- Per-chromosome inES files: ~25-50 MB each (compressed)
- Per-chromosome Rsb files: ~15-30 MB each (compressed)
- Genome-wide merged files: ~500 MB - 1 GB per population
- Peak extractions: ~5-20 MB total
- Summaries and percentiles: <1 MB

**Total storage:**
- GWAS peaks only: ~3-6 GB
- Full genome: ~8-12 GB

---

## Comparison with Other Methods

This EHH workflow complements other sweep detection methods:

- **SweepFinder2**: Frequency spectrum-based CLR statistic (best for hard sweeps)
- **SweeD**: Similar to SweepFinder2 but uses different SFS calculation
- **EHH (this workflow)**: Haplotype-based statistics (sensitive to recent sweeps, provides directionality with Rsb)

**Recommendation:** Run all three methods and look for **concordance**. Peaks showing high signals across multiple methods are most likely true selective sweeps.

---

## Questions?

See full documentation in:
- `code/12_Sweep_Detection/README.md` - Complete workflow overview
- `code/12_Sweep_Detection/ehh/methods.txt` - Methods description for papers (if created)
- `code/12_Sweep_Detection/ehh/ehh_merge_pheno.Rmd` - Existing peak-specific analysis (EHHS decay curves)
