# EHH & Selection: README

## What this module does (big picture)
This module takes a phased, biallelic SNP VCF and produces population-pair selection summaries and figures centered on GWAS peaks. The core logic is:

1. **Polarize the VCF** using an outgroup to infer the ancestral allele and add it as `INFO/AA`.  
2. **Split the polarized VCF** into per-population (or population-set) VCFs based on curated sample lists.  
3. **Optionally split** those VCFs by chromosome (or regions) to make downstream `rehh` operations fast and reproducible.  
4. **Load per-chromosome polarized VCFs in R (`rehh`)** and, for each GWAS peak and BP/TP population pair, compute:
   - **EHHS at the core SNP** (per population)  
   - **inES** (integrated EHH signals) around the peak  
   - **Rsb** (cross-population standardized iES difference)  
5. **Merge/plot** these results with figure panels and integrated summaries in `ehh_merge_pheno.Rmd`.

The result is a consistent, end-to-end path from the raw (but phased) variant data to interpretable selection signals around the GWAS peaks for the focal BP vs. TP comparisons.

---

## Repository layout (relevant bits)

```
code/
  12_Sweep_Detection/
    ehh/
      polarize_vcf.sh
      split_vcfs.sh
      split_vcf_by_chrom.sh
      analyze_peak_pair_integrated.R
      ehh_merge_pheno.Rmd

input_data/
  12_Sweep_Detection/
    *.txt
    split_up/

output_data/
  07_Phasing/
    ${phased_vcf}

  12_Sweep_Detection/
    AA.tsv.gz
    MIC4273_...polarized.vcf.gz
    lists/*.samples.txt
    vcfs/*.polarized.vcf.gz
    vcfs/by_chrom/*.vcf.gz
```

---

## Prerequisites
- **BCFtools** ≥ 1.13  
- **R** with packages: `rehh`, `tidyverse`, `patchwork`
- **Phased, biallelic SNP VCF** with consistent sample names
- **Outgroup samples** (default prefix: `PSZA_`)

---

## Quickstart

```bash
# 1) Polarize the full VCF using the outgroup
bash code/12_Sweep_Detection/ehh/polarize_vcf.sh

# 2) Split the polarized VCF into per-set VCFs
bash code/12_Sweep_Detection/ehh/split_vcfs.sh

# 3) (Optional) split by chromosome
bash code/12_Sweep_Detection/ehh/split_vcf_by_chrom.sh

# 4) Run the R notebook
R -e 'rmarkdown::render("code/12_Sweep_Detection/ehh/ehh_merge_pheno.Rmd")'
```

---

## Step-by-step

### 1. Polarize the VCF
Infers the ancestral allele from the outgroup and writes to `INFO/AA`.

### 2. Split by population sets
Creates VCFs per sample list (e.g., BP1_TP1).

### 3. Split by chromosome
Optionally subset per-chromosome for efficiency.

### 4. Analyze and plot in R
`ehh/ehh_merge_pheno.Rmd` loads per-chromosome polarized VCFs, runs `rehh`, and produces figures.

---

## Troubleshooting
- Ensure `INFO/AA` exists for polarization.
- Check sample name matches between lists and VCF.
- Edit chromosome targets in `split_vcf_by_chrom.sh`.
- Use per-chromosome subsets for large datasets.

---

---

## SweeD Genome-Wide Sweep Detection (Recommended)

In addition to the peak-specific EHH analysis above, this module includes a **SweeD** workflow for genome-wide sweep detection to establish background CLR distributions.

### Overview

**Purpose:** Calculate genome-wide Composite Likelihood Ratio (CLR) statistics to determine how "sweepy" GWAS peaks are relative to the rest of the genome.

**Method:** Site frequency spectrum (SFS)-based CLR approach (same statistical framework as SweepFinder2, but faster and more numerically stable)

**Populations:** BP1, TP1, BP2, TP2, BP3 (5 populations)

**Advantages over SweepFinder2:**
- **21x faster** sequential execution
- **Multi-threaded** parallel computation (8 threads per job)
- **Native VCF support** - no conversion step needed
- **Checkpoint/resume** capability for interrupted runs
- **Numerically more stable** floating-point operations

### Quick Start

See `code/12_Sweep_Detection/sweed/SWEED_QUICKSTART.md` for the complete guide.

**Basic workflow:**

```bash
# 1. Verify VCFs (optional, ~5 min)
bash code/12_Sweep_Detection/sweed/01_prepare_vcfs.sh

# 2. Calculate background SFS (~30 min)
bash code/12_Sweep_Detection/sweed/02_calculate_background_sfs.sh

# 3. Submit CLR scan jobs (~0.5-1 hour wall time)
bash code/12_Sweep_Detection/sweed/03_submit_sweed.sh

# 4. Merge results (~10 min)
bash code/12_Sweep_Detection/sweed/04_merge_clr_results.sh

# 5. Extract peak CLR (~5 min)
Rscript code/12_Sweep_Detection/sweed/05_extract_peak_clr.R

# 6. Generate report (~20 min)
bash code/12_Sweep_Detection/sweed/06_render_sweed_report.sh
```

### Expected Runtime

- **GWAS peaks (6 chromosomes):** 1.5-2.5 hours (60-70% faster than SweepFinder2)
- **Full genome (25 chromosomes):** 2-4 hours

### Prerequisites

- **SweeD container:** `images/sweed.sif` (build from `Dockerfile.sweed`)
- **Polarized VCFs:** Already created by `ehh/polarize_vcf.sh`
- **R** with tidyverse, knitr, rmarkdown packages

### Output Files

```
output_data/12_Sweep_Detection/sweed/
├── background/                 # Genome-wide empirical SFS
├── clr_results/               # Per-chromosome CLR results
├── merged/                    # Genome-wide CLR distributions
├── peak_clr_values.txt        # CLR within GWAS peaks
├── peak_clr_summary.txt       # Summary statistics
├── sweed_analysis.html        # Final HTML report
├── clr_manhattan.pdf          # Genome-wide plots
└── clr_zoomed_peaks.pdf       # Peak-specific plots
```

### When to Use SweeD vs SweepFinder2

**Use SweeD for:**
- New analyses (faster, more robust)
- Large datasets (>100k SNPs)
- Time-sensitive projects

**Use SweepFinder2 for:**
- Reproducing published results
- Direct comparison with older analyses
- Cross-validation

Both methods use the same CLR statistical framework and should give highly correlated results (r > 0.95).

---

## SweepFinder2 Genome-Wide Sweep Detection (Legacy)

This module also includes a **SweepFinder2** workflow for backward compatibility and validation.

### Overview

**Purpose:** Calculate genome-wide Composite Likelihood Ratio (CLR) statistics to determine how "sweepy" GWAS peaks are relative to the rest of the genome.

**Method:** Site frequency spectrum (SFS)-based CLR approach (does not require phasing)

**Populations:** BP1, TP1, BP2, TP2, BP3 (5 populations)

### Quick Start

See `code/12_Sweep_Detection/sweepfinder2/SWEEPFINDER2_QUICKSTART.md` for the complete guide.

### Prerequisites

- **SweepFinder2 container:** `images/sweepfinder2.sif` (build with `singularity build`)
- **Polarized VCFs:** Already created by `ehh/polarize_vcf.sh`
- **Python 3** with bcftools
- **R** with tidyverse, patchwork packages

### Workflow

```bash
# 0. Build SweepFinder2 container (one-time setup)
# See images/README.md for full instructions
#
# On local Mac:
cd images/
docker build --platform linux/amd64 -f Dockerfile.sweepfinder2 -t sweepfinder2 .
docker tag sweepfinder2 jasongallant/sweepfinder2  # Replace with your Docker Hub username
docker push jasongallant/sweepfinder2
#
# On HPCC:
cd images/
singularity build sweepfinder2.sif docker://jasongallant/sweepfinder2:latest

# 1. Convert polarized VCFs to SweepFinder2 format
bash code/12_Sweep_Detection/sweepfinder2/01_convert_vcf_to_sf2.sh
# Optional: --all-chromosomes flag to process all 25 chromosomes (default: GWAS peak chromosomes only)

# 2. Calculate genome-wide background SFS per population
bash code/12_Sweep_Detection/sweepfinder2/02_calculate_background_sfs.sh

# 3. Submit genome-wide CLR scan jobs (SLURM)
bash code/12_Sweep_Detection/sweepfinder2/03_submit_sweepfinder2.sh
# Optional: --all-chromosomes flag for full genome (default: GWAS peaks only)
# Monitor: squeue -u $USER

# 4. Merge per-chromosome CLR results (after jobs complete)
bash code/12_Sweep_Detection/sweepfinder2/04_merge_clr_results.sh

# 5. Extract GWAS peak CLR values
Rscript code/12_Sweep_Detection/sweepfinder2/05_extract_peak_clr.R

# 6. Generate analysis report
bash code/12_Sweep_Detection/sweepfinder2/06_render_sweepfinder2_report.sh
# Output: output_data/12_Sweep_Detection/sweepfinder2/sweepfinder2_analysis.html
```

### Computational Requirements

- **GWAS peaks only (6 chromosomes):** 30 jobs (5 populations × 6 chr), 2-3 hours wall time
- **Full genome (25 chromosomes):** 125 jobs (5 populations × 25 chr), 4-6 hours wall time
- **Memory:** 8 GB per job
- **Storage:** 3-6 GB total

### Output Files

```
output_data/12_Sweep_Detection/sweepfinder2/
├── sfs_input/              # SweepFinder2 format SFS files
├── background/             # Genome-wide background spectra
├── clr_results/            # Per-chromosome CLR results
├── merged/                 # Genome-wide CLR distributions
├── peak_clr_values.txt     # CLR within GWAS peaks
├── peak_clr_summary.txt    # Summary statistics
└── sweepfinder2_analysis.html  # Final report
```

### Key Outputs

1. **Genome-wide CLR distributions** - Empirical null model for each population
2. **Percentile rankings** - GWAS peak CLR values ranked against genome-wide background
3. **Statistical tests** - Wilcoxon, KS tests comparing peaks vs background
4. **Visualizations** - Manhattan plots, density plots, heatmaps

### Troubleshooting

- **Container not found:** Ensure `sweepfinder2_image` is set in `pkings_gwas.env` and container is built
- **No SFS files:** Check that `01_convert_vcf_to_sf2.sh` completed successfully
- **SLURM jobs failing:** Check logs in `output_data/slurm_logs/12_Sweep_Detection/`
- **Missing peak files:** Verify GWAS peaks BED file exists at `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed`

---

## Citation
- Gautier & Vitalis (2012, 2017), *rehh* R package for detecting selection signatures.
- Pavlidis et al. (2013), *SweeD* for genome-wide selective sweep detection (recommended).
- DeGiorgio et al. (2016), *SweepFinder2* for genome-wide selective sweep detection (legacy).

---

## Contact
Scripts authored and maintained by Jason Gallant, Michigan State University.
