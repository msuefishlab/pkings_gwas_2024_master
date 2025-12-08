# SweepFinder2 Quick Start Guide

## Overview

This workflow calculates genome-wide Composite Likelihood Ratio (CLR) statistics to determine how "sweepy" your GWAS peaks are relative to the rest of the genome.

**Populations:** BP1, TP1, BP2, TP2, BP3 (5 populations)
**GWAS Peaks:** 6 peaks (chr6, chr8, chr13, chr16, chr17, chr24)

---

## One-Time Setup: Build Container

### On your local Mac (M2):

```bash
cd images/
docker build --platform linux/amd64 -f Dockerfile.sweepfinder2 -t sweepfinder2 .
docker tag sweepfinder2 jasongallant/sweepfinder2
docker push jasongallant/sweepfinder2
```

### On HPCC:

```bash
cd images/
singularity build sweepfinder2.sif docker://jasongallant/sweepfinder2:latest
```

---

## Running the Analysis

### 1. Convert VCFs to SweepFinder2 Format (~1-2 hours)

```bash
source pkings_gwas.env
bash code/12_Sweep_Detection/01_convert_vcf_to_sf2.sh

# For full genome (optional):
# bash code/12_Sweep_Detection/01_convert_vcf_to_sf2.sh --all-chromosomes
```

**Output:** `output_data/12_Sweep_Detection/sweepfinder2/sfs_input/`

### 2. Calculate Background SFS (~30 minutes)

```bash
bash code/12_Sweep_Detection/02_calculate_background_sfs.sh
```

**Output:** `output_data/12_Sweep_Detection/sweepfinder2/background/`

### 3. Submit CLR Scan Jobs (~2-3 hours wall time)

```bash
bash code/12_Sweep_Detection/03_submit_sweepfinder2.sh

# For full genome (optional):
# bash code/12_Sweep_Detection/03_submit_sweepfinder2.sh --all-chromosomes

# Monitor jobs:
squeue -u $USER
```

**Jobs submitted:**
- GWAS peaks only: 30 jobs (5 populations × 6 chromosomes)
- Full genome: 125 jobs (5 populations × 25 chromosomes)

**Output:** `output_data/12_Sweep_Detection/sweepfinder2/clr_results/`

### 4. Merge Results (~10 minutes)

```bash
# After all jobs complete:
bash code/12_Sweep_Detection/04_merge_clr_results.sh
```

**Output:** `output_data/12_Sweep_Detection/sweepfinder2/merged/`

### 5. Extract Peak CLR Values (~5 minutes)

```bash
Rscript code/12_Sweep_Detection/05_extract_peak_clr.R
```

**Output:**
- `output_data/12_Sweep_Detection/sweepfinder2/peak_clr_values.txt`
- `output_data/12_Sweep_Detection/sweepfinder2/peak_clr_summary.txt`

### 6. Generate HTML Report (~20 minutes)

```bash
bash code/12_Sweep_Detection/06_render_sweepfinder2_report.sh
```

**Output:** `output_data/12_Sweep_Detection/sweepfinder2/sweepfinder2_analysis.html`

### 7. View Results

```bash
# Copy to your local machine or open in browser on HPCC
open output_data/12_Sweep_Detection/sweepfinder2/sweepfinder2_analysis.html
```

---

## What You'll Get

The HTML report includes:

1. **Genome-wide CLR distributions** - Percentiles (50th, 90th, 95th, 99th) per population
2. **Peak vs Background comparison** - Statistical tests (Wilcoxon, KS)
3. **Manhattan plots** - Genome-wide CLR with GWAS peaks highlighted
4. **Peak-specific plots** - CLR across each of 6 GWAS peaks
5. **BP vs TP comparison** - Population-specific selection patterns
6. **Percentile rankings** - How extreme each peak is
7. **Heatmaps** - Mean CLR by peak × population
8. **Summary statistics** - Key findings and interpretation

---

## Troubleshooting

**Container not found:**
```bash
# Check if container exists:
ls -lh images/sweepfinder2.sif

# If missing, rebuild from Docker Hub:
cd images/
singularity build sweepfinder2.sif docker://jasongallant/sweepfinder2:latest
```

**No SFS files after step 1:**
```bash
# Check conversion log:
cat output_data/12_Sweep_Detection/sweepfinder2/conversion.log

# Verify polarized VCFs exist:
ls output_data/12_Sweep_Detection/by_chrom/BP1.polarized.chr*.vcf.gz
```

**SLURM jobs failing:**
```bash
# Check job logs:
ls output_data/slurm_logs/12_Sweep_Detection/
tail output_data/slurm_logs/12_Sweep_Detection/sf2_BP1_chr6.log
```

**Missing peak BED file:**
```bash
# Verify GWAS peaks file exists:
ls output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed
```

---

## Total Time Estimate

- Container build (one-time): 15-30 minutes
- VCF conversion: 1-2 hours
- Background SFS: 30 minutes
- CLR scanning: 2-3 hours (GWAS peaks) or 4-6 hours (full genome)
- Merge + analysis + report: 35 minutes

**Total:** ~4-5 hours (GWAS peaks) or ~6-8 hours (full genome)

---

## File Sizes

- SFS input files: ~500 MB - 1 GB per population
- Background spectra: ~100 MB total
- CLR results: ~200-500 MB per population
- Merged files: ~500 MB - 1 GB per population
- **Total storage:** ~3-6 GB (GWAS peaks) or ~10-15 GB (full genome)

---

## Questions?

See full documentation in:
- `code/12_Sweep_Detection/README.md` - Complete workflow
- `code/12_Sweep_Detection/methods.txt` - Methods description for papers
- `images/README.md` - Container build instructions
