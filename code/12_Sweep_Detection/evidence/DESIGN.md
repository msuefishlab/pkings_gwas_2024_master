# Sweep evidence at the six GWAS peaks: pre-specified analysis design

**Status:** pre-specified. The commit that adds this file, with its hash and timestamp, is the record of what
was decided before the new statistics were computed. Changes after that commit must be made as further
commits, with a reason, and never by amending this one.

**Honest disclosure of what was already seen.** This is a re-analysis, not a first look. Before this design was
written we had already seen, for the six peaks:

- Rsb and inES percentiles against single SNPs and against equal-size genome windows
  (`output_data/12_Sweep_Detection/ehh/window_null/peak_window_null.tsv`, and the 2026-10-06 claims audit §1);
- SweepFinder2 CLR percentiles against same-size regions on the six peak chromosomes (exploratory, 2026-10-14);
- real iHS at Peaks 1 and 2 in BAM and APA (August 2025);
- the September 2025 π / Tajima's D / Fst mixed models.

So this is *not* a blind analysis, and we do not call it pre-registered. What the design fixes in advance is the
set of statistics, the search region, the window operator, the null, and the decision rule, so that these are
not chosen after seeing which combination favours a conclusion. The calibration in §6 is what makes the result
interpretable despite the prior looks: it measures how often the same rule fires on genome regions that are not
GWAS peaks.

## 1. Why this re-analysis

The published selection claims ("inES and iHS values at these BP-associated peaks fell above the 99th
percentile", "for all major BP-associated peaks, Rsb … exceeded the 99th percentile") rest on a comparison
between a window's *maximum* value and the distribution of *single SNPs*. A peak window holds 266–3,600 scored
SNPs, so its maximum is above the single-SNP 99th percentile almost regardless of whether a sweep occurred: in
fact 44–62% of non-peak 300-kb windows also contain at least one SNP above that threshold. The comparison
therefore does not test what the text claims it tests.

Three further problems are documented in the claims audit:

1. the statistic reported as "iHS" is inES standardized to mean 0 and SD 1 within each population, which is not
   iHS;
2. Fig. S7 plots population-level EHHS (`rehh::calc_ehhs`), not allele-specific EHH, so its legend's "BP
   (derived) versus TP (ancestral) allele" describes an analysis that was never run;
3. Fig. S7 shows BP3 at Peak 6, although the BP haplotype at Peak 6 is BP2-specific.

## 2. Units of analysis

Peaks are the extended GWAS peaks of decision D-002 (`input_data/10_Introgression/peaks.bed`): SNPs with
−log10 p > 3 merged within 75 kb. Index SNPs are the GWAS index SNPs in column 6 of that file.

| Peak | Gene | Chromosome | Interval | Index SNP | Lineage(s) in which the BP allele arose |
|---|---|---|---|---|---|
| 1 | c6orf132 | chr6 | 29,448,847–29,509,385 | 29,503,232 | BP1 |
| 2 | itih3 | chr8 | 15,771,736–15,805,670 | 15,786,956 | BP2, BP3 |
| 3 | cebpg/cebpa | chr13 | 24,628,382–24,929,106 | 24,780,908 | BP2, BP3 |
| 4 | cntnap5b | chr16 | 10,728,257–11,082,159 | 10,916,883 | BP1 |
| 5 | sptbn4b | chr17 | 11,494,062–11,781,225 | 11,692,887 | BP1 |
| 6 | naa40 | chr24 | 8,331,219–8,730,509 | 8,596,847 | BP2 |

Lineages and their matched references (`input_data/12_Sweep_Detection/split_up/`): BP1 = Bambomo (n = 42) with
TP1 = Ipassa + Bengue (n = 35); BP2 = Bikagala + Bavavela (n = 48) and BP3 = Biroundou (n = 16), both with
TP2 = Douvalou + Douengi (n = 31).

**Six tests, one per peak.** A peak's statistic is the maximum over its eligible lineages, and the null applies
the same maximum over the same lineages, so the multiple-lineage advantage is built into the null rather than
corrected for afterwards. BP2 and BP3 share BP alleles at Peaks 2 and 3, so they are not independent tests;
treating each peak as one test avoids claiming independence that does not exist. Per-lineage results are
reported in full alongside the six headline tests.

**Controls, reported but not part of the rule.** The matched TP population at each peak; the off-target BP
lineages; and BP1 at Peak 2, which the window null places at the 99.9th percentile although Peak 2's BP allele
did not arise in BP1 — a positive control for the machinery and a reminder that a high percentile alone does not
establish that a sweep drove the GWAS signal.

## 3. Applicability map, frozen before the statistics are computed

Step 1 (`01_index_snp_facts.R`) records, for every peak and every lineage, the frequency of the BP-associated
allele and whether it is the derived or the ancestral allele. The result is appended to this file and committed
*before* any statistic in §4 is computed, because it determines which statistics are even applicable:

- **iHS and nSL** compare haplotypes carrying the derived allele with those carrying the ancestral allele within
  one population. They lose power as the allele approaches fixation. They are declared applicable only where the
  BP allele frequency is ≤ 0.85 in that lineage.
- **Rsb, XP-EHH, SweepFinder2 CLR, π ratio and PBS** do not require intermediate frequency and apply everywhere.

Where a statistic is not applicable it is scored `NA` and is excluded from both the peak and the null, never
scored as a failure. A sweep that went to fixation must not be penalised for the fact that iHS cannot see it.

## 4. Statistics: one primary per evidence class

Three classes of evidence respond to different sweep ages and modes, and within each class the available
statistics are near-redundant (Rsb, XP-EHH, iHS and inES are all built from the same per-site integrated-EHH
quantities; SweeD and SweepFinder2 implement the same composite-likelihood ratio). Declaring one primary
statistic per class keeps the composite rule in §6 interpretable.

| Class | Primary statistic | Window operator | Rationale |
|---|---|---|---|
| H — haplotype structure | Rsb (BP versus matched TP) | mean over the best 20-kb window | standardized genome-wide, so it is not driven by the genome-wide BP/TP haplotype-length offset |
| S — site-frequency spectrum | SweepFinder2 CLR | maximum grid value in the region | detects sweeps after haplotypes have broken up; uses a genome-wide empirical background spectrum |
| D — diversity | π(BP)/π(TP) | minimum over the best 20-kb window | a sweep removes diversity in the swept lineage specifically |

**Secondary statistics**, reported in the evidence table and in the sensitivity file but not part of the
decision rule: XP-EHH; real iHS (`rehh::ihh2ihs`, freqbin 0.025); nSL (`scikit-allel`, standardized by
allele-frequency bin), which is less sensitive than iHS to variation in recombination rate; SweepFinder2 run on
the folded spectrum, as a check on polarization error; PBS for each BP lineage against its own TP and the other
TP; Tajima's D (descriptive only); and optionally SweeD on the same polarized input as SweepFinder2, as a
cross-implementation check rather than independent evidence.

**Excluded, with reasons.** Fay and Wu's H (redundant with the unfolded-spectrum CLR and maximally sensitive to
polarization error from a 12-individual outgroup); H12 and H2/H1 (distinguishing hard from soft sweeps is not
credible at n = 16 with unmeasured switch-error rates); allele age, haplotype networks, and
simulation-trained classifiers (out of scope for this revision); the SHAPEIT4-versus-WhatsHap switch-error check
(the per-sample WhatsHap outputs no longer exist on scratch).

## 5. Search region, window operator, and null

**Primary configuration, fixed in advance.** Search region: GWAS index SNP ± 100 kb. Window: the best
(highest for Rsb and CLR, lowest for the π ratio) 20-kb sub-window, stepping 5 kb, requiring at least 20 scored
SNPs. For the CLR the operator is instead the maximum grid value in the region, on the 10-kb grid.

The index SNP comes from the GWAS and is not chosen by any statistic computed here, so centring the search on it
does not use the sweep evidence twice.

**Operator identity.** The null applies the identical operator — same region size, same window size, same step,
same minimum SNP count, same maximum over lineages. The code asserts this rather than assuming it. Taking the
best of several window sizes or several region sizes would inflate the peak statistic relative to a null that
did not do the same; the alternatives in §7 are therefore reported as separate rows, each with its own matched
null, and never combined into the headline by taking the best.

**Null sets.** Pseudo-peaks are positions on all 25 chromosomes, excluding every peak ± 250 kb and assembly
gaps.

- **(a) Unmatched random pseudo-peaks**, at least 2,000. This answers the question of how many regions of this
  size anywhere in the genome look as sweep-like as the peaks do.
- **(b) Matched pseudo-peaks, the headline null.** Matched by decile on SNP density and on callable fraction
  (pixy's per-window count of sites with data). These two covariates affect every statistic here through the
  amount of information available in a window, and neither is a consequence of a sweep.

**Deliberately not matched on differentiation (user decision, 2026-10-14).** An earlier draft also matched on
the maximum BP–TP allele-frequency difference, on the grounds that the peaks were discovered *because* they are
differentiated. The objection that overrides it: a completed sweep *creates* high differentiation, so a null
matched on differentiation is built partly from other sweeps and from low-recombination regions where
differentiation accumulates by drift. Against such a null, "not supported" would mean only "no more sweep-like
than other highly differentiated regions", which is not the question being asked. Differentiation-matched
results are therefore reported as a sensitivity row (§7), where they place an upper bound on how much of the
signal could be explained by the ascertainment of the peaks.

## 6. Decision rule and calibration

For each peak, each primary statistic is converted to a percentile against the matched null (b).

- **Supported:** at least two of the three primary statistics at or above the 95th percentile.
- **Strong:** supported, and at least one primary statistic at or above the 99th percentile.
- **Composite score:** the second-highest of the three class percentiles — the value that the "two of three"
  rule turns on. The empirical p for a peak is the fraction of matched pseudo-peaks whose composite score is at
  least as high.

**Calibration is the headline result, not the percentiles.** The same rule is applied to every pseudo-peak in
both null sets, and we report the fraction that would be called "supported" and "strong". If the rule fires on
8% of random genome regions, then a peak that is merely "supported" is weak evidence, and the paper should say
so. This directly answers the earlier request to know how many outlier regions the same pipeline finds by
chance. We also report the Spearman correlation matrix among all statistics computed across the pseudo-peaks, so
that readers can judge how much independent information the "intersection of approaches" actually represents.

**Lineage specificity is a label, not a gate.** Each peak is labelled *BP-specific* (matched TP below the 95th
percentile), *shared* (TP also at or above the 95th), or *unsupported*. It is deliberately not a pass/fail
criterion: the TP population may carry its own sweep, and both populations rise together in regions of low
recombination, so requiring the TP to be low would discard real signals and, in a population of n = 31 or 35,
would mostly measure power.

## 7. Sensitivity rows (reported in full, never maxed into the headline)

Each is recomputed end to end, peak and null together: search region ± 50 kb; window size 40 kb; the D-002
extended peak interval in place of a fixed window; the unmatched null (a); the differentiation-matched null (c);
stratification of the null by a recombination proxy — the mean inES of the lineages *not* in the contrast, never
the matched TP, which may share the signal under test; SweepFinder2 on the folded spectrum; and every haplotype
statistic recomputed with all lineages down-sampled to n = 16, the size of BP3, with a fixed seed.

The down-sampled row needs a control row, added 2026-10-15 once the row was actually run. The only
down-sampled Rsb available is the one `10_genome_stats.R` recomputes from the saved `scan_hh` output, and it
standardizes Rsb across all 25 chromosomes at once, whereas the headline Rsb is standardized within each
chromosome by `ines2rsb` as each chromosome is scanned. Those two standardizations reorder regions relative to
one another between chromosomes, so comparing the down-sampled row directly with the headline would vary the
sample size and the standardization together and a changed call could not be attributed to either. The control
row (`_rsbgw`) is therefore the full sample with the genome-wide standardization, and the sample-size question
is answered by `_rsbgw` against `_n16`. This adds a comparison; it does not change the headline configuration
or the decision rule.

The down-sampling matters because the phasing (SHAPEIT4, no genetic map and no reference panel) has an unmeasured
switch-error rate that varies with sample size and heterozygosity, and switch errors shorten inferred
haplotypes. BP lineages carry genome-wide median inES 2–4× that of their TP references, which is exactly the
pattern a sample-size-dependent switch-error rate would produce. Any call that changes under down-sampling is
reported as such.

## 8. What the result can and cannot say

A peak that passes is consistent with a sweep on or near the BP-associated allele. It does not establish that
the sweep acted on the causal variant, nor that the sweep caused the phenotype. A peak that fails has not been
shown to be free of selection: a sweep that is old, soft, or from standing variation, or one in a region of high
recombination, may leave no detectable trace at this sample size. Both outcomes are reported with the
calibrated false-positive rate attached, and the manuscript text will follow the matrix rather than the
other way round.

## 9. Outputs

In `output_data/12_Sweep_Detection/evidence/`:

| File | Contents |
|---|---|
| `index_snp_facts.tsv` | BP-allele frequency and derived/ancestral state per peak per lineage; the applicability map |
| `evidence_matrix.tsv` | peak × statistic percentiles, composite score, empirical p, call, lineage-specificity label |
| `calibration.tsv` | false-positive rate of "supported" and "strong" in both null sets; per-statistic percentile uniformity; the Spearman correlation matrix |
| `sensitivity.tsv` | every row of §7, with the calls that change |

## 10. Provenance

Design written 2026-10-14. Scope and the choice not to match on differentiation were decided by J. Gallant on
2026-10-14; the demography-informed simulated null was deferred in the same discussion, with the reasoning
recorded in `private_data/decisions.md` (D-008). Analysis code lives in `code/12_Sweep_Detection/` on branch
`sweep_evidence`.
