# Drift test: fixed differences and the neutral expectation (D-009)

`DESIGN.md` is the pre-specified design and the authority for every choice here; this file only says how to
run it. The scripts in the parent directory (`prepare.sb`, `dadi_model.py`, `make_sfs.py`) are an older,
never-run single-population stub and are not part of this analysis.

## Run order (from the repository root on the HPCC)

| Step | Command | Produces |
|---|---|---|
| 0 | `sbatch --wrap` around `00_prepare_inputs.py --vcf_samples <bcftools query -l>` (gwas_tools.sif) | `input_data/11_Demographic_Modeling/{groups,group_sizes,pairs,chunks}.tsv`, `masks/{exons,exclude}.bed` |
| 1 | `sbatch --export=root=$(pwd) --array=0-102%40 code/11_Demographic_Modeling/fixed_differences/01_site_table.sb` | `site_table/{snps,invariant_hist,invariant_sums,index_snps}.<chunk>.*` |
| 2 | `python3 code/11_Demographic_Modeling/fixed_differences/02_check_site_table.py` | `checks/site_table_checks.tsv`; exits non-zero if the gate fails |

| 3–6 | `bash code/11_Demographic_Modeling/fixed_differences/run_to_spectra.sh` | spectra, projections, deserts, `sfs/<pair>.pkl`, `tail/index_units.tsv` |
| 7 | `run_fits.sh fit`, then `refine 2`, `refine 3` (each after the last finishes) | `fits/<pair>.<model>.data.seed*.tsv` |
| 7a | `run_fits.sh nested 1` after round 2, `nested 2` after round 3 (Addendum 2a; no-op if nothing is below its sub-model) | more seed files |
| 8 | `run_fits.sh collect <last_jobid>` | `fits/best_fits.tsv`, `fits/lrt.tsv` |
| 9 | `run_fits.sh tail <collect_jobid>` | `tail/{models,ladder,unit_results}.tsv` without bootstrap bounds |
| 10 | `run_fits.sh boot` (reads `tail/models.tsv` for the admissible best model), then `run_fits.sh tail <boot_jobid>` | final `tail/*.tsv` |

Set `FD_OUT=output_data/11_Demographic_Modeling/fixed_differences_agree` for the headline (Deviation 1) run; unset, the scripts write the *szaboi*-only run, now a sensitivity row. Report and figure follow DESIGN.md §9–§10.

## Things that will bite

- **Zero-coverage genotypes.** GATK writes `0/0` with DP = 0 for uncovered samples at invariant sites. `01_site_table.sb` sets them to missing before counting. Without that step, comparisons (and so L) run 2.6–5% high against pixy.
- **Outgroup.** It is the 11 P. szaboi of the popsfile. PSZA_455 is P. hopkinsi, even though its ID has the `PSZA_` prefix. Do not use `output_data/12_Sweep_Detection/AA.tsv.gz`: it came from the imputed phased VCF and included PSZA_455.
- **Two samplings of each creek.** `w*` groups are whole populations, as pixy used them; these are the headline. `p*` groups are the GWAS wobble-excluded lists; these are a sensitivity row. Both are computed in the same pass.
- **Input VCF.** It is on scratch (`${sweep_scratch}/03_QC/${allsites_vcf}`). The research-storage copies of this repository are unreadable; work in the rescue clone and push.
