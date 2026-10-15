# Sweep evidence: calibrated multi-method re-assessment of the six GWAS peaks

`DESIGN.md` is the pre-specified design and the authority for every choice here; this file only
says how to run it. Read DESIGN.md first — in particular §5 (the null) and §6 (the decision
rule and its calibration).

## The problem this replaces

The published selection claims compare a peak window's **maximum** inES or Rsb against the
distribution of **single SNPs**. A peak window holds 266–3,600 scored SNPs, so its maximum
clears the single-SNP 99th percentile almost regardless of whether a sweep occurred: 44–62% of
non-peak 300-kb windows contain such a SNP too. The comparison does not test what the text says
it tests.

Here every statistic is instead reduced by one fixed operator and compared against several
thousand non-peak regions of the same size, reduced by the identical operator, and the
decision rule is then applied to those same regions to measure how often it fires by chance.

## Run order

Everything runs inside `gwas_tools.sif` (bcftools, pixy, R with data.table) except the rehh and
SweepFinder2 steps, which use `rehh.sif` and `sweepfinder2.sif`. All paths are relative to the
repository root.

| Step | Command | Produces |
|---|---|---|
| 0 | `sbatch code/12_Sweep_Detection/00_fetch_inputs.sb` | containers, reference, the all-sites VCF (to scratch) |
| 1 | `sbatch code/12_Sweep_Detection/ehh/00_prepare_vcfs.sb` | polarized whole-genome VCF (or use the arrays below) |
| 1a | `sbatch code/12_Sweep_Detection/ehh/01a_split_sets_array.sb` | per-set whole-genome VCFs |
| 1b | `sbatch code/12_Sweep_Detection/ehh/01b_split_by_chrom_array.sb` | per-chromosome VCFs, 8 sets × 25 |
| 2 | `sbatch code/12_Sweep_Detection/ehh/14_ehh_scan_array.sb` | inES, Rsb, and the full `scan_hh` RDS |
| 3 | `Rscript code/12_Sweep_Detection/ehh/10_genome_stats.R` | real iHS, XP-EHH, Rsb check |
| 4 | `sbatch .../sweepfinder2/10_sf2_convert.sb` → `11_sf2_background.sb` → `12_sf2_scan.sb`, then `bash .../13_sf2_merge.sh` | CLR, all 25 chromosomes |
| 5 | `sbatch code/12_Sweep_Detection/pixy_lineages/run_pixy_lineages.sb` | π, dxy, Fst for the five lineage groups |
| 6 | `Rscript code/12_Sweep_Detection/evidence/01_index_snp_facts.R` | BP-allele frequencies, the applicability map |
| 7 | `Rscript .../evidence/02_region_stats.R` and `02b_region_callable.R` | per-region statistics and matching covariates |
| 8 | `Rscript .../evidence/03_evidence_matrix.R` | the evidence matrix and its calibration |
| 9 | `bash .../evidence/04_run_sensitivity.sh` | every sensitivity row of DESIGN.md §7 |

Steps 2–5 are independent of one another and can run at the same time; 1b must finish first.
Most scripts skip work that is already present, so a partly completed stage can be re-run.

## Files

- `region_stats_lib.R` — the search-region operator, the statistic loaders, and
  `assert_same_operator()`, which refuses to compare a peak and a null reduced differently.
  `build_regions()` defines the peaks and the pseudo-peak grid once, so every script scores
  exactly the same regions.
- `01_index_snp_facts.R` — BP-allele frequency and derived/ancestral state per peak per
  lineage; writes the applicability map that decides where iHS and nSL may be used.
- `02_region_stats.R` — applies the operator to the peaks and ~4,400 pseudo-peaks. Takes
  `--stats`, `--half_width`, `--window`, `--suffix`.
- `02b_region_callable.R` — callable sites (a matching covariate) and BP–TP Fst (used only by
  the differentiation sensitivity row).
- `03_evidence_matrix.R` — percentiles against the matched null, the two-of-three rule, the
  composite score, empirical p, the false-positive rate of the rule, and the correlation
  matrix among statistics. Takes `--suffix`, `--out_suffix`, `--null_mode`.
- `04_run_sensitivity.sh` — re-runs the reduction under each alternative configuration and
  prints the calls side by side.

## Reading the output

`evidence_matrix.tsv` carries, per peak and statistic, the value, its percentile against the
matched null, the same statistic in the matched TP population, and the call. **The percentile
on its own is not the result.** Read it next to `calibration.tsv`, which gives the fraction of
non-peak regions that the same rule calls "supported" and "strong". A peak that clears a bar
that 8% of random regions also clear is weak evidence, and the write-up should say so.

`statistic_correlations.tsv` shows how much the statistics agree across non-peak regions. Where
two statistics correlate strongly, their agreement at a peak is close to one observation, not
two.

## Things that will bite

- **Do not write into `output_data/12_Sweep_Detection/ehh/merged_published/`.** Those are the
  published scans restored from the S3 backup and are the reproduction target. New scans merge
  into `ehh/merged/`; the loader prefers that directory and falls back to the published one.
- The column earlier scripts wrote as `iHS` is inES standardized to mean 0 and SD 1, which is
  not iHS. It is now written as `inES_z`. Real iHS comes from `10_genome_stats.R`.
- `-lg` in SweepFinder2 is a spacing in base pairs, not a count. The published run's comment
  says 1 kb; its output is spaced at 10 kb.
- pixy must stay pinned at 1.2.10.beta2 (`images/Dockerfile`). Current master is 2.0.0.dev2 and
  fails on import. `pixy_lineages/validate_pixy.R` checks the rebuild against the published run
  through BP3, whose membership is identical to population BIR.
- Large intermediates live on scratch and are symlinked into `output_data/12_Sweep_Detection/`.
  They are gitignored; do not commit them.
