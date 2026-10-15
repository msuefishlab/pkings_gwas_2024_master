# Can drift explain the BP alleles at the GWAS peaks? Pre-specified design (D-009)

**Status:** pre-specified. The commit that adds this file, with its hash and timestamp (branch `drift_test`,
pushed to GitHub), records what was decided before any statistic in §4–§6 was computed. Changes after that
commit are made as further commits, each with a reason, and are listed under "Deviations" at the end.

**What had already been seen.** This design was written after the following were known:
- per-lineage BP-allele frequencies at the six index SNPs (`12_Sweep_Detection/evidence/index_snp_facts.tsv`);
- genome-wide pixy π, dxy, Fst and Tajima's D per population (`output_data/05_PopGen`);
- the D-008 sweep evidence matrix and calibration (only Peaks 3 and 6 in BP2 supported individually);
- quantities computed while planning:
  - da/π-based split times T/2N ≈ 0.14 (BAM–APA) and 0.27–0.45 (BIK, BAVA and BIR against DOV);
  - Fst ≈ 0.04 within the pooled TP references;
  - the exon-mask footprints (exons 11%, ±2 kb 49%, ±5 kb 63%, ±10 kb 76% of chr1–25);
- the group sizes, and the record counts, format and pixy-π agreement of two 1-Mb smoke tests (chr9:20–21 Mb, chr1:10–11 Mb; QC only, no fixed-difference counts).

No count of fixed differences, no joint SFS and no model fit had been looked at.

## 1. Question

Biphasic (BP) EODs evolved repeatedly in isolated *P. kingsleyae* populations. Earlier work disagreed on whether drift or selection explains this (drift in Gallant et al. 2011; selection in Picq et al. 2020). Genome-wide π shows large long-term Ne in every population, and the BP populations split from their TP neighbours recently relative to Ne. This analysis asks one question for each GWAS index SNP in each BP population:

> **How many SNPs as strongly BP-biased as this one would neutral drift have produced, anywhere in the genome, under the population history fitted from the genome itself?**

## 2. What this can and cannot show (written before any number)

- **Fixation alone proves nothing.** Every BAM fish is BP, so whatever causes BP must be near-fixed there whether drift or selection put it there. The test asks how often drift alone produces a difference this extreme, not whether one exists.
- **Standing variation is already in the null.** The fitted two-population history includes old polymorphisms that sort differently in the two daughter populations. So E (below) is the null for "drift sorted an old variant", and it covers the Peak 1 allele, which also occurs in BIR at 0.34.
  - The time argument (elapsed T/2Ne against the ~4Ne generations a new neutral mutation needs to fix) can only rule out "a new mutation drifted to fixation". It does not replace E.
- **The number of possible routes to BP (K) is unknown.** If K different sites could each make a fish BP, drift needs only one of them to reach this frequency, so its chance is about K × E / L. Here E is the expected number of such sites and L the number of callable sites. A broad mutational target (the Discussion; Reviewer 3's loss-of-function point) raises K and favours drift.
- **What we report about K.** For every unit we give K_needed = L / E: how many BP-capable sites the genome would need for drift to be expected to produce one. We also give it as a fraction of all callable sites. The reader can then judge whether a mutational target that large is plausible.
- **A K-independent check** (§6.4): drift with many routes predicts that the BP member of each pair is the lineage that drifted more.
- **What "selection" means here:** "not drift alone under the fitted class of histories". The test does not identify the selected trait or the target variant. "Drift readily makes this" is reported as "not evidence for selection", never as "neutral".

## 3. Units, samples and data

**Test units (11):** BAM × {P1, P4, P5}; BIK × {P2, P3, P6}; BAVA × {P2, P3, P6}; BIR × {P2, P3}. Each population is compared with its TP reference: BAM vs TP1 = APA + BEN; BIK, BAVA and BIR vs TP2 = DOV + DOG. BIK and BAVA are tested separately (Fst 0.24 between them).

Index SNPs come from column 6 of `input_data/10_Introgression/peaks.bed`: P1 chr6:29503232, P2 chr8:15786956, P3 chr13:24780908, P4 chr16:10916883, P5 chr17:11692887, P6 chr24:8596847.

**Headline sample: whole populations** (`input_data/05_PopGen/popsfile.txt`, wobble/mixed fish included), as pixy used them. This was decided by J. Gallant on 2026-10-15, reversing an earlier choice, for two reasons:
- Selecting fish by phenotype shifts allele frequencies only at the phenotype loci, which inflates the index SNP's differentiation relative to the genome; that bias works against drift.
- The demographic model describes populations.

The GWAS wobble-excluded lists are a sensitivity row. Group sizes: wBAM 67, wAPA 30, wBEN 17, wBIK 30, wBAVA 23, wBIR 16, wDOV 16, wDOG 15.

**Data:** the unphased all-sites VCF (`$allsites_vcf`; PASS sites only, invariant sites included), read in one `bcftools +fill-tags -S` pass (`01_site_table.sb`).
- **Ancestral allele:** the allele fixed in the 11 *P. szaboi* of the popsfile, with outgroup AN ≥ 8, otherwise unpolarized. PSZA_455 is *P. hopkinsi* and is excluded. *Superseded by Deviation 1 (outgroup agreement).*
- The older `AA.tsv.gz` is not used. It came from the SHAPEIT-imputed VCF, where missing outgroup calls were filled from mostly *P. kingsleyae* haplotypes, and it included PSZA_455.
- **Zero-coverage genotypes are missing.** GATK writes `0/0` with DP = 0 for uncovered samples at invariant sites. They are set to missing (`bcftools +setGT -i 'FMT/DP<1'`), which is also how pixy counts them; checked to reproduce pixy's comparison count exactly on one window.
- **Callable site for a pair:** AN ≥ the projection size in both populations and outgroup callable. SNPs are additionally required to be biallelic SNVs with a polarized ancestral allele. L (callable sites) uses the same rules.
- **Paralog filter:** SNPs with ExcHet < 0.001 in any whole population are dropped from data and L alike.

**Masks used for fitting:**
- Headline: exons excluded (all transcript types, from the RefSeq GFF), plus GWAS extended peaks ±250 kb and the chr4 inversion-like polymorphism (chr4:40.88–42.77 Mb). Only chr1–25 are used.
- Exons only, not a wider flank: a stricter mask fits a larger Ne, which predicts less drift near genes, and all six peaks sit near genes.
- Exons ±5 kb is a sensitivity row. Background selection is handled by the local-Ne correction E_hetNe (§4), not by a bigger mask.

**Call-desert mask** (fitting, L_mask and L_scan alike). The all-sites VCF keeps invariant sites fully called in some blocks where variant genotypes are mostly missing. An example is chr9:20–21 Mb: 987,288 invariant records and 908 SNP records, none of them in the 51M phased SNP set. These are probably blocks where invariant genotypes were written as zero-coverage `0/0` calls, or collapsed repeats where SNP genotypes failed a depth filter that invariant sites never saw. The DP = 0 rule should remove most of them, and the mask catches the rest.
- Such blocks add sites that the GWAS could not scan and that carry no SNPs.
- Pre-specified rule: per pair, a 100-kb bin is excluded when its callable SNPs per callable site fall below 10% of the genome-wide median for that pair.
- Excluded bins count as unscanned, so they are removed from L_scan too.
- For E_hetNe, bins whose outgroup divergence is below 25% of the genome median get f = 1.

**Projection:** hypergeometric, as in `dadi.Spectrum.from_data_dict`.
- Per population, n = the 10th percentile of AN over masked polarized SNPs, rounded down to an even number and capped at 60 haplotypes.
- Fallback, pre-specified: if one model evaluation takes more than 1 h, use n = 40.

## 4. Statistics

The statistics are per unit (BP population, TP reference, index SNP); i and j are derived-allele counts in BP and TP.

- **Δ_idx** = derived frequency in BP − derived frequency in TP at the index SNP, from raw counts in the headline sample. The BP allele must be derived (re-checked with the new ancestral calls).
- **Tail region R:** projected cells with i/n_BP − j/n_TP ≥ Δ_idx. For a fixed difference, R is the corner (n_BP, 0).
- **O**: observed SNPs in R. It is reported as O_mask (fitting mask), and as O_genome (all callable sites outside the peak's own ±250 kb).
- **E_mask**: the fitted model's expected count in R (θ-scaled, with the misidentification parameter). E_genome = E_mask × L_scan / L_mask, where L_scan = all callable sites on chr1–25.
  - P(drift makes at least one such site) ≤ E_genome by Markov's inequality, whatever the linkage. So "E ≪ 1" is a rigorous bound, not a heuristic.
- **E_hetNe** (background-selection-aware E):
  - Per 1-Mb window w, f_w = (π_off / d_OUT)_w divided by the same ratio genome-wide. π_off is the mean π of the four whole populations not in the pair; d_OUT is their mean divergence to the outgroup.
  - q is recomputed at local Ne × f (dadi units: T/f, m·f, θ·f), interpolated on a grid of f from 0.2 to 1.5; then E_hetNe = Σ_w L_w · q(f_w).
  - Also reported: f_peak, the peak's own window value, and f\*, the local Ne at which the peak region alone would be expected to produce one such site.
- **τ_BP, τ_TP**: each branch's total drift since the split, ∫dt / 2N(t) from the fitted model. F = 1 − e^(−τ), and T/2N_anc.
- **K_needed** = L_scan / E_genome (§2).
- **Ladder:** O and E_mask at Δ ∈ {0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0}, outside the exclusions.
- **Asymmetry:** observed and model-expected counts of BP-derived (Δ ≥ 0.95) and TP-derived (Δ ≤ −0.95) differences.
- **Parallel count**, for P2, P3 and P6 (reported, not part of the rule): genome-wide SNPs where the same derived allele is at least as frequent as the index SNP in each of BIK, BAVA (and BIR for P2 and P3), and at most as frequent in TP2.
- **LD-proxy row**, for P4 and P5 (not fixed): E for the least-differentiated SNP with r² ≥ 0.8 to the index SNP in the BP population, from the phased VCF. Δ then no longer depends on which SNP in the block is causal.

## 5. Models

- dadi 2.4.3 (`images/dadi_tools.sif`). Every model is wrapped with the ancestral-misidentification parameter. Grid `pts = [n+10, n+20, n+30]` with extrapolation.
- Models; BP is population 1, so in M3 s is the BP founder fraction:
  - **M0** split, no migration;
  - **M1** split, symmetric migration (`split_mig`);
  - **M2** split, asymmetric migration (`split_asym_mig`);
  - **M3** `IM`: founder fraction s, exponential growth, asymmetric migration.
- **Fitting:** 20 perturbed starts per model × pair, using `Inference.opt`.
  - Converged = the best three runs lie within 1 log-likelihood unit, with parameters within 5%.
  - Uncertainty comes from Godambe (GIM) with 100 bootstrap spectra from 1-Mb blocks. 5-Mb blocks are a sensitivity row.
  - Godambe-adjusted LRTs (M0 ⊂ M1 ⊂ M2) are reported for information. No AIC: dadi's likelihood is composite.
- **Admissible model:** converged, misidentification rate ≤ 0.05, and it passes the adequacy gate: at the most extreme ladder threshold with E ≥ 10, O/E is within [1/3, 3].
- **Best model:** the admissible model with the highest log-likelihood. 100 bootstrap refits under the best model and under M3 give the confidence intervals for E, τ and T/2N.
- Full model set on the 4 headline pairs and the 2 TP–TP control pairs (wAPA–wBEN, wDOV–wDOG). Sensitivity configurations get M3 plus the best model.

## 6. Decision rule

### 6.1 E_robust (drift-generous by construction)

E_robust is the largest of:
- the upper 95% bootstrap bound of E_genome under the best model;
- the point estimate of E_genome under every admissible model;
- E_hetNe under the best model.

### 6.2 Tiers

| E_robust | Reading |
|---|---|
| ≤ 0.005 | **Drift implausible**, even allowing for all 11 units. |
| 0.005–0.05 | **Drift implausible for this unit** (P ≤ 0.05 for this unit alone). |
| 0.05–1 | **Unusual under drift.** Suggestive only. |
| 1–20 | **Drift makes a few such sites.** Differentiation alone does not show selection. |
| > 20 | **Drift readily makes this.** Not evidence for selection, and not evidence against it. |

### 6.3 Overrides

- **Model fails at the extremes:** if E_robust ≤ 1 but O_genome > max(5, 3·E_genome), the genome holds more such sites than the model allows. The empirical count O becomes the headline for that unit, and the tier is read from O, which is drift-generous. If the TP–TP controls show the same excess, the cause is model inadequacy, not BP-specific selection.
- **No admissible model for a pair:** that pair's units are "model-dependent". Report O and the empirical rank only.
- **Model-dependent units:** a unit whose tier differs between admissible models is labelled "model-dependent"; its tier is the most drift-generous one.
- **Pooled versus single TP references:** if a call changes by more than one tier between the pooled TP reference and a single TP population (APA, BEN, DOV or DOG), the more drift-generous call is the headline.

### 6.4 Drift-exposure check (K-independent; reported next to the tiers)

- Drift with many routes to BP predicts that the BP population of each pair drifted more than its TP reference: τ_BP > τ_TP, and more BP-derived than TP-derived near-fixed differences beyond the model's expectation.
- **Pre-specified reading:**
  - If τ_BP ≤ τ_TP in at least 3 of 4 headline pairs, drift-with-many-routes cannot explain why these populations became BP and their neighbours did not. The text says so.
  - If τ_BP > τ_TP in at least 3 of 4, the drift explanation stays open regardless of the tiers, and the text says that instead.

### 6.5 Controls

- **Off-target index SNPs** are reported for every pair, for example P1 in BIR–TP2 (0.34 vs 0), and P2, P3 and P6 in BAM–TP1. They are expected in the "drift readily" tier or to be unremarkable.
- **TP–TP pairs** are expected to have O/E ≈ 1 at the corners.

## 7. Sensitivity rows (reported in full; never maxed into the headline except under §6.3)

1. GWAS wobble-excluded lists (pBAM … pTP2).
2. Single TP references: APA, BEN, DOV, DOG.
3. BIK + BAVA pooled (pBP2/wBP2), for continuity with D-008.
4. Mask exons ±5 kb.
5. Folded SFS; its tail E_dir + E_reverse is an upper bound on the directional E.
6. 5-Mb bootstrap blocks.
7. Singletons and doubletons masked in the fit.
8. Paralog filter off.
9. LD-proxy SNP for P4 and P5.

## 8. Phase B: GWAS-aware null for the D-008 sweep statistics (conditional)

Phase B runs only if the D-008 lineage pairs have at least 200 non-peak, 75-kb-merged clusters of BP-derived near-fixed SNPs (≥ 0.95 / ≤ 0.05). The pairs use pure lineages, as D-008 does: pBAM/pTP1, pBP2/pTP2 and pBIR/pTP2.
- If a lineage falls short, relax to 0.9 / 0.1. If it still falls short, Phase B is not run for that lineage. That is reported as a result: drift-made fixed differences are too rare to form a null.
- **Centres:** each cluster's most differentiated SNP, excluding peaks ±250 kb.
- **Code:** the D-008 operator is reused unchanged through a new `--null_mode fixed_diff`. It is added to the D-008 DESIGN.md as a dated row, never as an amendment.
- **Reading:** a peak supported here is more sweep-like than other drift-sorted fixed differences. The null contains other sweeps, so it is conservative.

## 9. Gates and verification

1. **Site table** (`02_check_site_table.py`), which must pass before §4:
   - whole-population π per population and chromosome within 0.5% of pixy, and dxy for the six check pairs within 0.5%;
   - SNP records compared with the 51,073,851-site phased SNP set. The all-sites VCF also carries SNPs that the phasing inputs dropped (e.g. chr9:20–21 Mb), so its count is expected to be higher; the overlap and the excess are reported, not gated;
   - index-SNP frequencies (pure lineage groups) within 0.05 of `index_snp_facts.tsv`;
   - the BP allele derived at all six index SNPs.
2. **Spectra:**
   - folding the unfolded spectrum equals the directly folded spectrum;
   - totals equal the retained sites;
   - the count-table path matches `from_data_dict` on a 1-Mb subset;
   - block sums equal the genome-wide spectrum.
3. **Models:**
   - convergence as in §5;
   - model-implied π, dxy and Fst match the observed values;
   - the sign of Tajima's D matches for BIK (+1.57) and BAM;
   - Poisson parameter recovery (20 replicates) recovers parameters and E within their intervals;
   - corner cells change by less than 5% between pts and about 2×pts. If not, E is reported as "numerically unreliable".
4. msprime is not used for the headline. An msprime check of the corner cells and cluster counts is added only if an E_robust falls between 0.05 and 20, as a recorded deviation.

## 10. Outputs (`output_data/11_Demographic_Modeling/fixed_differences/`)

- `site_table/`: per-chunk SNP tables, invariant histograms and sums.
- `checks/site_table_checks.tsv`.
- `fixed_diff_counts.tsv`, `fixed_diff_clusters.<pair>.bed`, `ladder.tsv`.
- `sfs/`: data and bootstrap spectra, and `projections.tsv`.
- `fits/model_fits.tsv`, with residual plots.
- `tail/unit_results.tsv`: per unit, Δ_idx, O, E (all variants), E_robust, tier, K_needed, f_peak, f\*, τ and the controls.
- `RESULTS.md` and the report HTML, with the plain-language row per unit and the tail figure: genome-wide count of sites at least as extreme against Δ, observed vs model, with the index SNPs marked.

**Sentence template for each unit:**

> "Since [population] split from its TP neighbours, drift has had time to move a typical 50% allele by about ±[√(F·0.25)]. The elapsed time is [T/2N] × 2Ne generations, [x]% of the time a new neutral mutation needs to fix, so any neutral fixed difference must come from variation that already existed. Under the best-fitting history, drift would create [E] (95% CI) sites in the genome as strongly BP-biased as this SNP. We observe [O] outside the GWAS peaks. Verdict: [tier]."

## 11. Provenance

- Design written 2026-10-15, after an independent design review. Scope (this revision), whole-population headline, and Phase B were decided by J. Gallant on 2026-10-15. Decision record: `private_data/decisions.md` D-009.
- Code: `code/11_Demographic_Modeling/fixed_differences/` on branch `drift_test`.

## Deviations

### Deviation 1 (2026-10-15): ancestral allele by agreement of six outer-clade species

**What was seen first.** The first set of fits used the *P. szaboi*-only rule (§3). It returned an ancestral-misidentification rate of 3–10% in every pair, including both TP–TP controls (APA–BEN about 10%). §5 makes such models inadmissible above 5%. No E value, tier or unit result had been computed when this deviation was decided. The observed fixed-difference counts and the parallel counts had been seen; neither depends on the labelling.

**Why.** The species tree (`output_data/04_Phylogenomics/...cf.tree`) places *P. szaboi* inside an outer clade with *P. magnostipes*, *P. parsegnis*, *P. hopkinsi*, *P. gabonensis* and *P. vadamans*. Gene concordance at the deep nodes is low (about 12–40%), which points to extensive incomplete lineage sorting. Many *P. kingsleyae* polymorphisms are therefore old, and a single outgroup species fixed for one allele may simply have sorted it by chance. Species sort independently, so requiring several to agree makes such errors much less likely. Reference bias in the outgroup calls may also contribute; agreement among species does not remove that, which is why the misidentification rate is re-estimated (see the readings below).

**New rule.** The six outer-clade species are separate outgroup groups (OUT_PSZA 11, OUT_PMAG 20, OUT_PPAR 7, OUT_PHOP 5, OUT_PGAB 4, OUT_PVAD 2 fish; PSZA_455 is *P. hopkinsi* and counts under PHOP).
- A SNP is polarized only if at least 2 species have data (AN ≥ 2 each), they carry at least 8 alleles in total, and every species with data is fixed for the same allele.
- An invariant site is outgroup-callable under the same data rule.
- The curvifrons group (*P. curvifrons*, *P. ntotom*, *P. longicaudatus*, *P.* sp. SN2) is not used. It is closer to *P. kingsleyae* and shares excess derived alleles with it (D-006).

**Unchanged.** Everything else in this design: masks, projections, models, the 5% admissibility rule, tiers, sensitivity rows.

**Kept for comparison.** The *P. szaboi*-only run is kept in `output_data/11_Demographic_Modeling/fixed_differences/` and reported as a sensitivity row. The new run writes to `output_data/11_Demographic_Modeling/fixed_differences_agree/` (`FD_OUT`).

**How the re-estimate will be read.** *Superseded by Addendum 1a, before any rerun fit existed.*

Decided by J. Gallant, 2026-10-15.

#### Addendum 1a (2026-10-15, before any fit of the rerun): reading the rate, and the folded headline

**Correction.** The reading above, that a rate staying above 5% means the models are too simple, is wrong. Two labelling errors survive outgroup agreement:
- **Mapping bias.** All six species are mapped to the same *P. kingsleyae* reference and sit about equally far from it, so they share any bias toward the reference allele.
- **Shared sorting.** The six species share a common ancestor after their split from *P. kingsleyae*. An old polymorphism that sorted in that ancestor makes all six agree on the wrong allele.

So a rate that stays high could be either leftover mislabelling or too-simple models. These fits alone cannot tell the two apart.

**Evidence already in hand** (`output_data/03_QC/ref_bias/ref_bias_per_sample.tsv`, D-005). The share of reference reads at heterozygous calls is 0.543 in *P. kingsleyae* and 0.544 in *P. szaboi*; in *P. magnostipes*, *parsegnis*, *hopkinsi* and *gabonensis* it is 0.566–0.582. Outgroups have more missing calls (16–17.5% vs 12.6%).
- *P. szaboi* shows no more reference skew than the ingroup, so mapping bias is unlikely to explain the *szaboi*-only rate. The other outgroups show a slight extra skew.
- Mapping bias mostly produces missing calls rather than wrong ones: a fish carrying only the non-reference allele produces no reference reads.

**Rule** (J. Gallant, 2026-10-15). The label-free (folded) analysis is built for every pair from every SNP, polarized or not. It is fitted with M0f–M3f: the same demographies without a misidentification parameter. Its tail counts SNPs at least as differentiated as the index SNP *in either direction*, so its E can only exceed the directional E, the drift-generous direction.
- If at least one unfolded model of a pair is admissible (§5, including misidentification ≤ 5%), the unfolded analysis is the headline for that pair, as designed.
- If no unfolded model is admissible, the folded analysis is the headline, subject to the same convergence and adequacy gates. The unfolded results are reported alongside.
- The rate is reported with both explanations named. A rate that falls to about 1–2% is read as mislabelling fixed. A rate that stays high is read as "mislabelling or model inadequacy, not distinguished", and the folded headline does not depend on which.
