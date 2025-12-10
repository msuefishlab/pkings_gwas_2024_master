# EHH Analysis Pipeline

## Overview

Extended Haplotype Homozygosity (EHH) analysis for paired populations to detect selection signatures and compute genome-wide percentile rankings for GWAS peaks.

**Paired populations:**
- BP1_TP1
- BP2_TP2
- BP3_TP2

**Statistics computed:**
- **inES** - Integrated Extended Haplotype Statistic (raw EHH measure per population)
- **iHS** - Integrated Haplotype Score (standardized inES, mean=0, SD=1)
- **Rsb** - Relative Selection Between populations (cross-population comparison)

**Key feature:** This pipeline uses a **paired population approach** that guarantees marker alignment between BP and TP populations, matching the proven peak analysis workflow.

---

## Prerequisites

- Phased VCF from step 07_Phasing
- Sample lists in `input_data/12_Sweep_Detection/split_up/` (BP1.txt, BP2.txt, BP3.txt, TP1.txt, TP2.txt)
- Outgroup samples (PSZA_*) for ancestral allele inference

---

## Workflow

### Step 1: Prepare VCFs (one-time setup)

#### 1a. Add ancestral allele annotation

```bash
bash code/12_Sweep_Detection/ehh/01_polarize_vcf.sh
```

**What it does:**
- Uses outgroup (PSZA samples) to infer ancestral allele
- Adds INFO/AA field to VCF
- Required for polarized EHH analysis

**Output:** `MIC4273_...polarized.vcf.gz`

**Time:** ~30-60 minutes

---

#### 1b. Create paired population VCFs

```bash
bash code/12_Sweep_Detection/ehh/02_split_by_sample_sets.sh
```

**What it does:**
- Reads sample lists from `input_data/12_Sweep_Detection/split_up/`
- Creates paired VCFs (BP1_TP1, BP2_TP2, BP3_TP2)
- Each VCF contains both populations for that comparison

**Output:**
- `output_data/12_Sweep_Detection/vcfs/BP1_TP1.polarized.vcf.gz`
- `output_data/12_Sweep_Detection/vcfs/BP2_TP2.polarized.vcf.gz`
- `output_data/12_Sweep_Detection/vcfs/BP3_TP2.polarized.vcf.gz`

**Time:** ~1-2 hours

---

#### 1c. Split by chromosome

```bash
bash code/12_Sweep_Detection/ehh/03_split_by_chromosome.sh

# For GWAS peak chromosomes only (default):
bash code/12_Sweep_Detection/ehh/03_split_by_chromosome.sh

# OR for all 25 chromosomes:
bash code/12_Sweep_Detection/ehh/03_split_by_chromosome.sh --all-chromosomes
```

**What it does:**
- Splits each paired VCF into per-chromosome files
- Default: GWAS peak chromosomes (chr6, chr8, chr13, chr16, chr17, chr24)
- Optional: All 25 chromosomes

**Output:**
- `output_data/12_Sweep_Detection/by_chrom/BP1_TP1.polarized.chr6.vcf.gz`
- `output_data/12_Sweep_Detection/by_chrom/BP1_TP1.polarized.chr8.vcf.gz`
- ... (18 VCF files for GWAS peaks, or 75 for all chromosomes)

**Time:** ~10-30 minutes

---

### Step 2: Run genome-wide EHH scans

```bash
# Submit integrated EHH scan jobs (GWAS peaks)
bash code/12_Sweep_Detection/ehh/04_submit_ehh_scans.sh

# Monitor jobs
squeue -u $USER

# For all 25 chromosomes (optional):
# bash code/12_Sweep_Detection/ehh/04_submit_ehh_scans.sh --all-chromosomes
```

**What it does:**
For each pair × chromosome, computes in a **single integrated pass**:
1. Load paired VCF
2. Apply MAF filter to combined dataset (ensures marker alignment)
3. Split into BP and TP subsets on same marker set
4. Compute inES/iHS for BP population
5. Compute inES/iHS for TP population
6. Compute Rsb between populations

**Jobs submitted:**
- GWAS peaks: **18 jobs** (3 pairs × 6 chromosomes)
- All chromosomes: **75 jobs** (3 pairs × 25 chromosomes)

**Resources per job:**
- Time: 3 hours
- Memory: 16 GB
- CPUs: 1

**Output (3 files per job):**
- `output_data/12_Sweep_Detection/ehh/scans/BP1_TP1.BP.chr6.ines.txt.gz` (BP1 inES/iHS)
- `output_data/12_Sweep_Detection/ehh/scans/BP1_TP1.TP.chr6.ines.txt.gz` (TP1 inES/iHS)
- `output_data/12_Sweep_Detection/ehh/scans/BP1_TP1.chr6.rsb.txt.gz` (Rsb)

**Total files:**
- GWAS peaks: 18 jobs × 3 files = **54 files**
- All chromosomes: 75 jobs × 3 files = **225 files**

**Expected wall time:** ~3 hours (if sufficient nodes available)

---

### Step 3: Merge results

```bash
# After all jobs complete successfully
bash code/12_Sweep_Detection/ehh/05_merge_ehh_results.sh
```

**What it does:**
- Merges per-chromosome inES files into genome-wide files per population
- Merges per-chromosome Rsb files into genome-wide files per comparison
- Deduplicates TP2 (appears in both BP2_TP2 and BP3_TP2)

**Output:**
```
output_data/12_Sweep_Detection/ehh/merged/
├── BP1.genome_wide_ines.txt.gz     # BP1 population inES/iHS
├── BP2.genome_wide_ines.txt.gz     # BP2 population inES/iHS
├── BP3.genome_wide_ines.txt.gz     # BP3 population inES/iHS
├── TP1.genome_wide_ines.txt.gz     # TP1 population inES/iHS
├── TP2.genome_wide_ines.txt.gz     # TP2 population inES/iHS (deduplicated)
├── BP1_TP1.genome_wide_rsb.txt.gz  # BP1 vs TP1 Rsb
├── BP2_TP2.genome_wide_rsb.txt.gz  # BP2 vs TP2 Rsb
├── BP3_TP2.genome_wide_rsb.txt.gz  # BP3 vs TP2 Rsb
└── merge.log                        # Detailed merge log
```

**Total:** 5 inES files + 3 Rsb files = **8 merged genome-wide files**

**Time:** ~10-30 minutes

---

### Step 4: Extract peak statistics

```bash
# Source environment to get access to rehh_image variable
source pkings_gwas.env

# Run with singularity
singularity exec --bind ${root}:/project_root ${rehh_image} \
  Rscript /project_root/code/12_Sweep_Detection/ehh/06_extract_peak_ehh.R
```

**What it does:**
1. Calculate genome-wide summary statistics (mean, SD, percentiles)
2. Extract EHH values within GWAS peaks
3. Calculate percentile rankings and empirical p-values

**Output:**
- `genome_wide_ines_summary.txt` - Genome-wide percentiles for inES/iHS
- `genome_wide_rsb_summary.txt` - Genome-wide percentiles for Rsb
- `peak_ines_values.txt` - All SNPs within GWAS peaks (inES/iHS)
- `peak_ines_summary.txt` - Summary per peak (mean, median, max, SD)
- `peak_rsb_values.txt` - All SNPs within GWAS peaks (Rsb)
- `peak_rsb_summary.txt` - Summary per peak
- `peak_ines_percentiles.txt` - Percentile rank of each peak's max value
- `peak_rsb_percentiles.txt` - Percentile rank for Rsb
- `peak_ines_empirical_pvalues.txt` - Statistical significance
- `peak_rsb_empirical_pvalues.txt` - Statistical significance

**Time:** ~15 minutes

**Example interpretation:**
> "Peak 2 in BP1 has max |iHS| = 2.345, which is in the 99.87th percentile genome-wide (only 0.13% of SNPs have higher |iHS|). Empirical p-value = 0.0013 (genome-wide significant)."

---

### Step 5: Generate report

```bash
bash code/12_Sweep_Detection/ehh/07_render_ehh_report.sh
```

**What it does:**
- Renders R Markdown analysis to HTML
- Creates visualizations (Manhattan plots, heatmaps, distributions)

**Output:**
- `output_data/12_Sweep_Detection/ehh/ehh_analysis.html` (main report)

**Report includes:**
1. Genome-wide summary statistics - Mean, SD, percentiles for all populations/comparisons
2. Percentile heatmaps - Visual summary of peak rankings (rows=peaks, columns=populations/comparisons)
3. Empirical p-values - Statistical significance table
4. Distribution comparison plots - Genome-wide density (gray) with peak values overlaid (red rug)
5. Manhattan plots - Genome-wide landscape with GWAS peaks highlighted (blue), percentile thresholds (dashed lines)
6. Peak summary tables - Top-ranked peaks sorted by percentile

**Generated figures (saved as PDFs):**
- `ihs_percentile_heatmap.pdf` - Heatmap of iHS percentiles
- `rsb_percentile_heatmap.pdf` - Heatmap of Rsb percentiles
- `ihs_distribution_comparison.pdf` - Genome-wide vs peak distributions (iHS)
- `rsb_distribution_comparison.pdf` - Genome-wide vs peak distributions (Rsb)
- `ihs_manhattan.pdf` - Genome-wide Manhattan plot (iHS)
- `rsb_manhattan.pdf` - Genome-wide Manhattan plot (Rsb)

**Time:** ~20 minutes

**View results:**
```bash
# Open in browser
open output_data/12_Sweep_Detection/ehh/ehh_analysis.html
```

---

## File Structure

```
code/12_Sweep_Detection/ehh/
├── 01_polarize_vcf.sh              # Add ancestral allele
├── 02_split_by_sample_sets.sh      # Create paired VCFs
├── 03_split_by_chromosome.sh       # Split into per-chr VCFs
├── 04_submit_ehh_scans.sh          # Submit EHH scan jobs
├── 04_run_ehh_scan.sb              # SLURM batch script
├── 05_merge_ehh_results.sh         # Merge per-chr results
├── 06_extract_peak_ehh.R           # Extract peak stats
├── 07_render_ehh_report.sh         # Generate HTML report
├── scan_ehh_integrated.R           # Core R scanner (called by 04_run_ehh_scan.sb)
└── analyze_peak_pair_integrated.R  # Helper functions
```

**Deprecated scripts** (kept for reference):
- `DEPRECATED_submit_ines_scans.sh`
- `DEPRECATED_submit_rsb_scans.sh`
- `DEPRECATED_run_ines_scan.sb`
- `DEPRECATED_run_rsb_scan.sb`
- `DEPRECATED_scan_ines_genome_wide.R`
- `DEPRECATED_scan_rsb_genome_wide.R`

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
- **High iHS + High Rsb** → Strong population-specific sweep at EOD locus
- **High iHS + Low Rsb** → Shared selection across populations
- **Low iHS + High Rsb** → Old sweep with population-specific allele frequency differences

---

## Troubleshooting

### Missing VCF files

```
ERROR: VCF file not found: /path/to/BP1_TP1.polarized.chr6.vcf.gz
```

**Solution:** Run `03_split_by_chromosome.sh` first to create per-chromosome VCFs.

```bash
bash code/12_Sweep_Detection/ehh/03_split_by_chromosome.sh
```

---

### Marker count mismatch

```
ERROR: Marker count mismatch after split!
```

**Solution:** This should **never happen** with the integrated pipeline because the `split_hh_by_groups()` function ensures both populations share the same marker set. If you see this error, it indicates a bug - please report it.

---

### Job failure due to memory

```
slurmstepd: error: Exceeded job memory limit
```

**Solution:** The integrated scan uses 16GB memory. If still failing:

1. Check which chromosome failed (large chromosomes like chr1, chr2 need more memory)
2. Increase memory in `04_run_ehh_scan.sb`:
   ```bash
   #SBATCH --mem=24GB  # Increase from 16GB
   ```
3. Re-submit failed jobs manually

---

### Job failure due to time limit

```
slurmstepd: error: *** JOB CANCELLED AT ... DUE TO TIME LIMIT ***
```

**Solution:** Large chromosomes (chr1, chr2) may need >3 hours. Increase time limit in `04_run_ehh_scan.sb`:

```bash
#SBATCH --time=04:00:00  # Increase from 03:00:00
```

Then re-submit failed jobs.

---

### Not all merged files created

```
WARNING: Not all merged files were created (6 / 8)
```

**Solution:** Check that all scan jobs completed successfully:

```bash
# Check merge log for details
cat output_data/12_Sweep_Detection/ehh/merged/merge.log

# Check for failed jobs in SLURM logs
ls output_data/slurm_logs/12_Sweep_Detection/ehh/
tail output_data/slurm_logs/12_Sweep_Detection/ehh/ehh_BP1_TP1_chr6.*.log

# Verify all chromosomes completed
ls output_data/12_Sweep_Detection/ehh/scans/ | grep BP1_TP1 | wc -l
# Should be 18 (6 chromosomes × 3 file types) for GWAS peaks

# Identify which chromosomes/pairs failed
ls -l output_data/12_Sweep_Detection/ehh/scans/ | grep -E "BP|TP|rsb"

# Re-submit failed jobs manually using sbatch
```

---

### Check job status

```bash
# Monitor running jobs
squeue -u $USER

# Watch jobs in real-time
watch -n 60 'squeue -u $USER'

# Cancel all EHH jobs (if needed)
scancel -u $USER -n ehh
```

---

### Verify prerequisite files

```bash
# Check that polarized VCF exists
ls output_data/12_Sweep_Detection/*.polarized.vcf.gz

# Check that paired VCFs exist
ls output_data/12_Sweep_Detection/vcfs/*.polarized.vcf.gz

# Check that per-chromosome VCFs exist
ls output_data/12_Sweep_Detection/by_chrom/*.polarized.chr*.vcf.gz
```

---

## Performance

### Time Estimates

**VCF Preparation (one-time):**
- 01_polarize_vcf.sh: ~30-60 minutes
- 02_split_by_sample_sets.sh: ~1-2 hours
- 03_split_by_chromosome.sh: ~10-30 minutes

**Analysis (GWAS peaks):**
- 04_submit_ehh_scans.sh: ~3 hours wall time (parallelized)
- 05_merge_ehh_results.sh: ~10-30 minutes
- 06_extract_peak_ehh.R: ~15 minutes
- 07_render_ehh_report.sh: ~20 minutes

**Total:** ~4-6 hours for GWAS peak chromosomes

**Total (all chromosomes):** ~6-10 hours for all 25 chromosomes

### Resource Usage

**Compared to old pipeline:**
- **62.5% fewer jobs** (48 → 18 for GWAS peaks)
- **50% less I/O** (loads each paired VCF once instead of twice)
- **55% fewer node-hours** (~100 → ~45 for GWAS peaks)
- **Guaranteed marker alignment** between BP and TP populations

### Disk Space

**Per-file estimates:**
- Per-chromosome inES files: ~25-50 MB each (compressed)
- Per-chromosome Rsb files: ~15-30 MB each (compressed)
- Genome-wide merged inES files: ~500 MB - 1 GB per population
- Genome-wide merged Rsb files: ~400 MB - 800 MB per comparison
- Peak extractions: ~5-20 MB total
- Summary files and percentiles: <1 MB

**Total storage:**
- GWAS peaks only: ~3-6 GB
- Full genome (all 25 chromosomes): ~8-12 GB

---

## Key Implementation Details

### Paired Population Approach

The pipeline processes populations in **pairs** (BP1_TP1, BP2_TP2, BP3_TP2) to ensure marker alignment:

1. Load paired VCF containing both populations
2. Apply MAF filter to **combined dataset** (before splitting)
3. Split into BP/TP subsets using `subset(hh_all, select.hap = indices)`
4. This guarantees **identical SNP positions** for both populations

**Why this matters:** If populations were filtered separately, different SNPs might be removed due to MAF differences, causing position mismatches when computing Rsb.

### Direct Rsb Computation

The pipeline computes `ines2rsb()` directly from `scan_hh()` result objects, not from saved files. This ensures:
- Exact position matching (no need for alignment logic)
- Preserves all metadata from scan objects
- Matches the proven peak analysis approach

### TP2 Deduplication

TP2 is intentionally used in both BP2_TP2 and BP3_TP2 comparisons (biological design). The merge script:
1. Processes TP2 from BP2_TP2 pair first
2. Appends TP2 from BP3_TP2 pair
3. Deduplicates by CHR:POSITION (keeping first occurrence)
4. Creates single `TP2.genome_wide_ines.txt.gz` file

---

## Comparison with Other Methods

This EHH workflow complements other sweep detection methods:

- **SweepFinder2**: Frequency spectrum-based CLR statistic (best for hard sweeps)
- **SweeD**: Similar to SweepFinder2 but uses different SFS calculation
- **EHH (this workflow)**: Haplotype-based statistics (sensitive to recent sweeps, provides directionality with Rsb)

**Recommendation:** Run all three methods and look for **concordance**. Peaks showing high signals across multiple methods are most likely true selective sweeps.

---

## Questions or Issues?

Contact: Jason Gallant Lab

**Common questions:**

- **Q: Can I process a subset of chromosomes?**
  - A: Yes, modify the `CHROMOSOMES` array in `04_submit_ehh_scans.sh` before running.

- **Q: What if I only want to re-run one failed chromosome?**
  - A: Use `sbatch` to manually submit a single job with the appropriate environment variables. See `04_submit_ehh_scans.sh` for the exact command format.

- **Q: How do I validate results match the peak analysis?**
  - A: Compare Rsb values for peak regions between `ehh_merge_pheno.Rmd` output and genome-wide merged files. They should be highly correlated (R² > 0.99).
