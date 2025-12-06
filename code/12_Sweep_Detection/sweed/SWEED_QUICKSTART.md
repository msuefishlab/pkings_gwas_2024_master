# SweeD Quick Start Guide

## Overview

SweeD is a fast, efficient tool for genome-wide selective sweep detection using Composite Likelihood Ratio (CLR) statistics. This guide walks through the complete SweeD workflow for the P. kingsleyae GWAS project.

**Key advantages over SweepFinder2:**
- **21x faster** sequential execution
- **Multi-threaded** parallel computation (8 threads per job)
- **Native VCF support** - no conversion step needed
- **Checkpoint/resume** capability for interrupted runs
- **Numerically more stable** floating-point operations

**Expected runtime:**
- GWAS peaks only (6 chromosomes): 1.5-2.5 hours
- Full genome (25 chromosomes): 2-4 hours

---

## Prerequisites

- **Polarized VCFs** with INFO/AA ancestral allele annotations
  - Created by `code/12_Sweep_Detection/polarize_vcf.sh`
  - Per-population VCFs from `split_vcfs.sh`
  - Optional: per-chromosome VCFs from `split_vcf_by_chrom.sh`

- **GWAS peaks BED file**
  - `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed`

- **SweeD container** (built from `images/Dockerfile.sweed`)

- **R packages**: tidyverse, knitr, rmarkdown

---

## One-Time Setup: Build SweeD Container

### On Local Mac (M2/ARM)

```bash
cd images/

# Build for linux/amd64 (HPCC architecture)
docker build --platform linux/amd64 -f Dockerfile.sweed -t sweed .

# Tag with your Docker Hub username
docker tag sweed <your-dockerhub-username>/sweed:latest

# Push to Docker Hub
docker push <your-dockerhub-username>/sweed:latest
```

### On HPCC (Linux)

```bash
cd images/

# Pull and convert to Singularity
singularity build sweed.sif docker://<your-dockerhub-username>/sweed:latest

# Verify
singularity exec sweed.sif SweeD -h
```

### Update Environment

Ensure `pkings_gwas.env` contains:
```bash
export sweed_image=${root}/images/sweed.sif
```

---

## Workflow: 6-Step Pipeline

### Step 1: Verify VCFs (Optional, ~5 minutes)

Validates that polarized VCFs exist with INFO/AA annotations:

```bash
# Source environment
source pkings_gwas.env

# Verify VCFs (GWAS peaks only)
bash code/12_Sweep_Detection/sweed/01_prepare_vcfs.sh

# OR verify all chromosomes
bash code/12_Sweep_Detection/sweed/01_prepare_vcfs.sh --all-chromosomes
```

**Output:**
- `output_data/12_Sweep_Detection/sweed/vcfs/vcf_manifest.txt` - List of VCFs to process
- `output_data/12_Sweep_Detection/sweed/vcf_validation.log` - Validation report

**What it checks:**
- ✅ VCF files exist for each population × chromosome
- ✅ VCFs have INFO/AA field (ancestral allele)
- ⚠️  Warns if AA annotations are missing

**Note:** Unlike SweepFinder2, SweeD reads VCF directly - **no conversion step needed!**

---

### Step 2: Calculate Background SFS (~30 minutes)

Generates genome-wide empirical site frequency spectrum (SFS) per population:

```bash
bash code/12_Sweep_Detection/sweed/02_calculate_background_sfs.sh
```

**What it does:**
- Runs SweeD genome-wide per population with `-osfs` flag
- Outputs empirical SFS for use as background null model
- Uses 8 threads for parallel computation

**Output:**
- `output_data/12_Sweep_Detection/sweed/background/BP1.genome_sfs.txt`
- `output_data/12_Sweep_Detection/sweed/background/TP1.genome_sfs.txt`
- `output_data/12_Sweep_Detection/sweed/background/BP2.genome_sfs.txt`
- `output_data/12_Sweep_Detection/sweed/background/TP2.genome_sfs.txt`
- `output_data/12_Sweep_Detection/sweed/background/BP3.genome_sfs.txt`

**Alternative:** Skip this step to use SweeD's theoretical SFS (constant population size) - faster but less realistic.

---

### Step 3: Submit CLR Scan Jobs (0.5-1 hour wall time)

Submits SLURM jobs for genome-wide CLR calculation:

```bash
# GWAS peaks only (30 jobs: 5 pops × 6 chromosomes)
bash code/12_Sweep_Detection/sweed/03_submit_sweed.sh

# OR full genome (125 jobs: 5 pops × 25 chromosomes)
bash code/12_Sweep_Detection/sweed/03_submit_sweed.sh --all-chromosomes
```

**Job resources:**
- CPUs: 8 per job (multi-threading)
- Memory: 16 GB per job
- Time limit: 4 hours
- Checkpoint: Every hour (can resume if interrupted)

**Monitor jobs:**
```bash
# Check status
squeue -u $USER

# Count running jobs
squeue -u $USER -n sweed | wc -l

# Watch job progress
watch -n 30 'squeue -u $USER -n sweed | wc -l'
```

**Cancel jobs if needed:**
```bash
# Cancel all SweeD jobs
scancel -u $USER -n sweed
```

**Output:**
- `output_data/12_Sweep_Detection/sweed/clr_results/BP1.chr6.sweed.txt`
- `output_data/12_Sweep_Detection/sweed/clr_results/BP1.chr8.sweed.txt`
- ... (one file per population × chromosome)
- Log files: `output_data/slurm_logs/12_Sweep_Detection/sweed_*.log`

**SweeD output format:**
```
// Position    Likelihood  Alpha
1000          2.456       0.0015
2000          1.234       0.0020
```

---

### Step 4: Merge CLR Results (~10 minutes)

Combines per-chromosome CLR files into genome-wide distributions:

```bash
bash code/12_Sweep_Detection/sweed/04_merge_clr_results.sh
```

**What it does:**
- Filters out SweeD comment lines (`//`)
- Adds chromosome column
- Merges all chromosomes per population
- Renames `Likelihood` → `CLR` for consistency

**Output:**
- `output_data/12_Sweep_Detection/sweed/merged/BP1.genome_wide_clr.txt`
- `output_data/12_Sweep_Detection/sweed/merged/TP1.genome_wide_clr.txt`
- `output_data/12_Sweep_Detection/sweed/merged/BP2.genome_wide_clr.txt`
- `output_data/12_Sweep_Detection/sweed/merged/TP2.genome_wide_clr.txt`
- `output_data/12_Sweep_Detection/sweed/merged/BP3.genome_wide_clr.txt`

**Merged format:**
```
chromosome  position  CLR    alpha
chr6        1000      2.456  0.0015
chr6        2000      1.234  0.0020
```

---

### Step 5: Extract GWAS Peak CLR (~5 minutes)

Extracts CLR values within the 6 GWAS peak regions and calculates summary statistics:

```bash
Rscript code/12_Sweep_Detection/sweed/05_extract_peak_clr.R
```

**What it does:**
- Loads genome-wide CLR distributions
- Filters to positions within GWAS peaks
- Calculates mean, median, max CLR per peak per population
- Identifies position of maximum CLR signal

**Output:**
- `output_data/12_Sweep_Detection/sweed/peak_clr_values.txt` - All CLR values within peaks
- `output_data/12_Sweep_Detection/sweed/peak_clr_summary.txt` - Summary statistics

**Peak summary columns:**
- `population`: BP1, TP1, BP2, TP2, BP3
- `peak`: Peak number (1-6)
- `chromosome`: chr6, chr8, chr13, chr16, chr17, chr24
- `n_positions`: Number of grid points in peak
- `mean_CLR`: Average CLR within peak
- `median_CLR`: Median CLR within peak
- `max_CLR`: Maximum CLR within peak
- `sd_CLR`: Standard deviation of CLR
- `position_max`: Position of maximum CLR

---

### Step 6: Generate HTML Report (~20 minutes)

Renders comprehensive analysis report with visualizations:

```bash
bash code/12_Sweep_Detection/sweed/06_render_sweed_report.sh
```

**Output:**
- `output_data/12_Sweep_Detection/sweed/sweed_analysis.html` - Interactive HTML report
- `output_data/12_Sweep_Detection/sweed/clr_manhattan.pdf` - Genome-wide Manhattan plot (PDF)
- `output_data/12_Sweep_Detection/sweed/clr_manhattan.png` - Genome-wide Manhattan plot (PNG)
- `output_data/12_Sweep_Detection/sweed/clr_zoomed_peaks.pdf` - Zoomed peak plots (PDF)
- `output_data/12_Sweep_Detection/sweed/clr_zoomed_peaks.png` - Zoomed peak plots (PNG)

**View report:**
```bash
# On local machine
open output_data/12_Sweep_Detection/sweed/sweed_analysis.html

# Or download from HPCC
scp hpcc:~/path/to/output_data/12_Sweep_Detection/sweed/sweed_analysis.html .
```

**Report sections:**
1. **Genome-wide CLR Manhattan plot** - All chromosomes, all populations
2. **Zoomed peak regions** - ±250kb windows around each GWAS peak
3. **Summary statistics** - Peak vs background CLR comparisons
4. **Per-peak CLR summary table** - Top 20 peaks by maximum CLR
5. **Percentile distributions** - CLR percentiles (50th, 90th, 95th, 99th)

---

## Troubleshooting

### Container Issues

**Problem:** SweeD container not found
```
ERROR: SweeD container not found: /path/to/images/sweed.sif
```

**Solution:**
```bash
# Check if container exists
ls -lh images/sweed.sif

# If missing, build it (see "One-Time Setup" above)
```

---

### VCF Issues

**Problem:** VCFs missing AA annotations
```
WARNING: X VCFs are missing ancestral allele (INFO/AA) annotations!
```

**Solution:**
```bash
# Run polarization workflow
bash code/12_Sweep_Detection/polarize_vcf.sh
bash code/12_Sweep_Detection/split_vcfs.sh
bash code/12_Sweep_Detection/split_vcf_by_chrom.sh
```

**Verify AA field:**
```bash
bcftools view -h output_data/12_Sweep_Detection/vcfs/BP1.polarized.vcf.gz | grep "INFO=<ID=AA"
```

---

### Job Failures

**Problem:** SweeD jobs failing

**Check logs:**
```bash
# View latest log
ls -lt output_data/slurm_logs/12_Sweep_Detection/sweed_*.log | head -1 | xargs cat

# Search for errors
grep -i error output_data/slurm_logs/12_Sweep_Detection/sweed_*.log
```

**Common issues:**
1. **Out of memory** - Increase `--mem` in `03_run_sweed.sb` from 16GB to 32GB
2. **VCF not found** - Verify VCF paths in Step 1
3. **Background SFS missing** - Re-run Step 2 or remove `-isfs` flag in `03_run_sweed.sb`

---

### Merge Errors

**Problem:** Merged files have wrong format

**Check SweeD output:**
```bash
# View raw SweeD output
head -20 output_data/12_Sweep_Detection/sweed/clr_results/BP1.chr6.sweed.txt
```

**Expected format:**
```
// Position    Likelihood  Alpha
// (comment lines)
1000          2.456       0.0015
```

**Solution:** Merge script filters `//` lines correctly - check grep command in `04_merge_clr_results.sh`

---

### R Errors

**Problem:** R script fails to load data

**Check file existence:**
```bash
# Verify merged files
ls -lh output_data/12_Sweep_Detection/sweed/merged/

# Verify peak files
ls -lh output_data/12_Sweep_Detection/sweed/peak_*.txt
```

**Check R packages:**
```r
# In R
library(tidyverse)
library(rprojroot)
```

---

## Comparing SweeD vs SweepFinder2

### CLR Correlation

After running both methods, compare CLR distributions:

```r
# In R
library(tidyverse)

# Load SweeD results
sweed_clr <- read_tsv("output_data/12_Sweep_Detection/sweed/merged/BP1.genome_wide_clr.txt")

# Load SweepFinder2 results
sf2_clr <- read_tsv("output_data/12_Sweep_Detection/sweepfinder2/merged/BP1.genome_wide_clr.txt")

# Match positions (allowing some tolerance)
merged <- inner_join(sweed_clr, sf2_clr, by = c("chromosome", "position"),
                     suffix = c("_sweed", "_sf2"))

# Calculate correlation
cor(merged$CLR_sweed, merged$CLR_sf2, use = "complete.obs")
# Expected: r > 0.95
```

### Runtime Comparison

**SweepFinder2:**
- Step 1 (VCF conversion): 1-2 hours
- Step 2 (Background SFS): 30 min
- Step 3 (CLR scan): 2-3 hours
- Steps 4-6: 35 min
- **Total: 4-5 hours**

**SweeD:**
- Step 1 (VCF validation): 5 min (NO CONVERSION!)
- Step 2 (Background SFS): 30 min
- Step 3 (CLR scan): 0.5-1 hour
- Steps 4-6: 35 min
- **Total: 1.5-2.5 hours**

**Speedup: 60-70% reduction**

---

## Performance Benchmarks

### GWAS Peaks (6 chromosomes, 5 populations)

- **Jobs:** 30 total
- **Runtime:** 1.5-2.5 hours
- **Storage:** ~2 GB
- **CPU-hours:** ~120 (8 CPUs × 0.5 hr × 30 jobs)

### Full Genome (25 chromosomes, 5 populations)

- **Jobs:** 125 total
- **Runtime:** 2-4 hours
- **Storage:** ~8 GB
- **CPU-hours:** ~500 (8 CPUs × 0.5 hr × 125 jobs)

---

## Next Steps

After completing the SweeD workflow:

1. **Review HTML report** - Check genome-wide CLR patterns and peak enrichment

2. **Identify "sweepy" peaks** - Peaks with CLR > 95th percentile genome-wide

3. **Integrate with haplotype analysis** - Compare CLR results with EHH/Rsb statistics from `ehh_merge_pheno.Rmd`

4. **Validate results** - Compare with SweepFinder2 outputs to ensure concordance

5. **Publication** - Use methods text from `sweed/methods.txt` and figures from report

---

## Citation

If you use SweeD in your analysis, please cite:

> Pavlidis P, Zivkovic D, Stamatakis A, Alachiotis N (2013) SweeD: Likelihood-Based Detection of Selective Sweeps in Thousands of Genomes. Molecular Biology and Evolution 30(9): 2224-2234. doi:10.1093/molbev/mst112

---

## Contact

For questions or issues with this workflow:
- **GitHub Issues:** https://github.com/alachins/sweed/issues
- **Lab Contact:** Jason Gallant Lab, Michigan State University
