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

Later steps (spectra, fits, tail expectations, report) are added as they are written. Each one follows DESIGN.md §4–§9.

## Things that will bite

- **Zero-coverage genotypes.** GATK writes `0/0` with DP = 0 for uncovered samples at invariant sites. `01_site_table.sb` sets them to missing before counting. Without that step, comparisons (and so L) run 2.6–5% high against pixy.
- **Outgroup.** It is the 11 P. szaboi of the popsfile. PSZA_455 is P. hopkinsi, even though its ID has the `PSZA_` prefix. Do not use `output_data/12_Sweep_Detection/AA.tsv.gz`: it came from the imputed phased VCF and included PSZA_455.
- **Two samplings of each creek.** `w*` groups are whole populations, as pixy used them; these are the headline. `p*` groups are the GWAS wobble-excluded lists; these are a sensitivity row. Both are computed in the same pass.
- **Input VCF.** It is on scratch (`${sweep_scratch}/03_QC/${allsites_vcf}`). The research-storage copies of this repository are unreadable; work in the rescue clone and push.
