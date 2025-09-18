## Purpose
In this module, we ask "Do GWAS peaks in BP1, BP2, and BP3 show localized excess allele sharing consistent with introgression?"  

### Whole-genome D-statistics with Dsuite (Dtrios)
`01_submit_dsuite.sh` submits chromosome-by-chromosome Dsuite Dtrios jobs inside a Singularity image, using:

	•   a rooted species/population tree (your pops.tree)
	•	a SETS.txt mapping samples → populations
	•	SNP-only VCFs sliced by chromosome (via bcftools view)

 Dsuite computes (per trio P1,P2,P3 + outgroup from tree):

	•	D (Patterson’s D; ABBA–BABA test): excess of ABBA vs BABA site patterns.
	•	Sign: D > 0 → excess sharing between P3 and P2 (relative to P1); D < 0 → excess sharing between P3 and P1.
	•	Z and p: block-jackknife significance tests for D.
	•	f4-ratio: rough admixture proportion proxy (how much of P3’s ancestry might trace to P2 vs P1 given the topology).
	•	BBAA/ABBA/BABA counts: diagnostics for how many informative sites contributed.

Interpretation: Under ILS alone, ABBA≈BABA and D≈0; strong |D| suggests introgression (or model violations).

`DsuiteAnalysis.Rmd` reads the BBAA/Dtrios outputs, standardizes column names, and produces summaries/plots (e.g., ranking trios by Z, exploring directions of gene flow, comparing f4-ratio across trios).

•	Is introgression is present at all, and between which population pairs?

•	Is introgression widespread (many trios, many chromosomes) or sparse/focal (a few trios, weak genome-wide signal).

### Windowing to candidate peaks and re-testing (Dtrios on peaks)

`02_make_peak_window.sh` takes a BED of GWAS peaks and pads each peak by a user-specified half-window (±K kb) to create an analysis window BED.

`03_submit_dsuite_peaks.sh` then restricts VCFs to those windows (via bcftools with -R/regions) and re-runs Dsuite Dtrios only in those intervals.


•	Do D, Z, and f4-ratio increase in peak windows vs. genome-wide baseline?

•	Are signals consistent across multiple peak regions or specific to one/few?

•	If introgression is enriched at candidate loci, consistent with the idea that a phenotype-associated allele might have moved between populations.

•	If peak windows do not show elevated D/f4 relative to background (despite strong association), that argues against introgression as the source of the shared signal there and is more consistent with parallel/de-novo changes (or balancing selection within lineages, etc.).

### Localizing signals with sliding windows (Dinvestigate)

`04_submit_dinvestigate.sh` + `run_dsuite_investigate.sb` launch Dsuite Dinvestigate for one or more specific trios Dinvestigate slides along the genome (or your regions) in fixed windows.

What Dinvestigate computes per window:

	•	D again (local), plus f_d, f_dM, and d_f (window-wise normalizations of f-like estimators; commonly used to pinpoint introgressed tracts).
	•	f_d: tends to be conservative and is bounded (useful for highlighting peaks of introgression).
	•	f_dM: a modified version aiming to reduce biases across allele frequency spectra and topologies.
	•	d_f: another normalization; you’re using it as a complementary lens—peaks that agree across f_d and f_dM are especially compelling.

`DsuiteInvestigateAnalysis.Rmd` ingests the windowed outputs, joins them to  peak windows (labels), and plots tracks of D / f_d / f_dM / d_f across coordinates.

	• Is the signal sharp and co-localized to your candidate interval (suggestive of introgressed haplotypes at/near the causal locus), or diffuse (suggestive of broader demographic admixture or structure)?

	•Consistency across statistics: Peaks supported by D + f_d + f_dM are stronger evidence than peaks supported by only one metric.

### Glossary

• **D (ABBA–BABA)** - Excess shared derived alleles. D≈0 under ILS; |D|>0 suggests introgression (or topology/ascertainment issues). Sign indicates which sister shares with P3 more (positive → P2–P3; negative → P1–P3 given your ordering).

• **Z (for D)** - Significance via block-jackknife. Typical informal thresholds are |Z|≥3, but you contextualize this with multiple testing and consistency across chromosomes.

• **f4-ratio (Dtrios)** - A rough admixture proportion estimate for the focal trio’s direction (bounded between 0 and 1 in well-behaved cases). Use it to compare magnitude of putative gene flow among trios.

• **f_d, f_dM, d_f (Dinvestigate)** - Window-level analogs/normalizations designed to pinpoint local introgressed tracts; peaks co-localized with your GWAS hits are especially interesting.

•	**BBAA / ABBA / BABA counts** - Numbers of sites informing each pattern; low counts can produce unstable estimates—important for QC.

