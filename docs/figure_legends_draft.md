# Figure legends — DRAFT (discovery set: GSE67528, GSE124308, GSE142670)

Numbers are from `logs/07b_figures_revised.log` and `logs/08_main_figures.log` (24 Nov 2026).
Items in [brackets] still need checking or a decision.

---

## Figure 1. Study flow and dataset overview
[Not yet made: PRISMA-style flow of the search and screening, plus the dataset table (cohort, model, stage, donors per group, samples). Built after the formal search.]

---

## Figure 2. Direction of ASD-associated expression changes replicates across independent cohorts.
Leave-one-cohort-out replication. For each held-out cohort, the 500 genes with the strongest combined signal in the remaining cohort(s) (weighted Stouffer z, weights √donors) were selected, and the share changing in the same direction in the held-out cohort was recorded. Grey, models adjusted for donor only (`dream`, donor as random effect); green, additionally adjusted for a per-sample neuronal marker score. Error bars, exact 95% binomial confidence intervals; dashed line, 50% expected by chance. Fewer than 500 genes were testable where the held-out cohort did not express all selected genes.
**(A) Neurons.** Unadjusted agreement was 62.9% (n = 498 genes; one-sided binomial p = 5.3 × 10⁻⁹), 61.5% (n = 465; p = 4.0 × 10⁻⁷) and 61.3% (n = 452; p = 9.2 × 10⁻⁷) for GSE67528, GSE124308 and GSE142670, respectively; after neuron-score adjustment, 68.7% (p = 1.8 × 10⁻¹⁷), 63.2% (p = 4.9 × 10⁻⁹) and 57.9% (p = 4.8 × 10⁻⁴). ‡Sensitivity analysis replacing GSE142670 mixed-culture neurons with FACS-sorted neurons from the same study: agreement was 32.6% (n = 500), significantly *below* chance (two-sided p = 5.4 × 10⁻¹⁵), i.e. a reversed direction [discussed in Results; cause unresolved].
**(B) Neural progenitor cells.** Two cohorts only; each cohort's top 500 genes were tested in the other. Agreement when GSE67528 was held out: 62.6% (n = 423; p = 1.1 × 10⁻⁷); when GSE142670 was held out: 55.6% (n = 360; p = 0.020). Neuron-score adjustment gave 62.6% (p = 6.4 × 10⁻⁸) and 56.6% (p = 7.4 × 10⁻³).

---

## Figure 3. Pathway themes across cohorts and stages.
Gene-set enrichment of selected Hallmark, GO Biological Process and Reactome sets. Genes were ranked by the moderated t statistic of the ASD-vs-control coefficient in each cohort, or by the weighted Stouffer z for the pooled columns (shaded; not independent of the cohort columns). Fill, fgsea normalised enrichment score (NES; red, higher in ASD; blue, lower); dot size, fgsea −log₁₀ p (capped at 10). fgsea p values do not account for inter-gene correlation; black outlines mark sets that were also significant by the correlation-adjusted competitive test cameraPR (inter-gene correlation 0.05; Benjamini–Hochberg FDR < 0.05 across the 12 sets within each column). Cilia sets were lower in ASD in GSE67528 and GSE124308 neurons, in GSE142670 NPCs and in both pooled rankings, but not in GSE142670 neurons. MYC targets V1 were higher in ASD NPCs in both cohorts (cameraPR FDR 4.7 × 10⁻³ and 3.9 × 10⁻⁷). †Interferon sets were pre-specified because of the earlier analysis of these data; the pooled-neuron interferon-α result (FDR 0.024) is driven by GSE124308 alone, the cohort in which iPSC source is confounded with diagnosis (Fig. S1B).

---

## Figure 4. The cilia signal is partly explained by FOXJ1 expression.
cameraPR test of the curated motile-cilia set (top; the pre-specified primary test for H2) and of the curated primary-cilia set (bottom), before and after adding a covariate to the donor-aware model. Bars show signed −log₁₀ p (negative, lower in ASD); the shaded band marks p > 0.05. One-sided p values for the motile-cilia set being lower in ASD (unadjusted → + FOXJ1): GSE67528 neurons 2.1 × 10⁻⁹ → 6.0 × 10⁻⁴; GSE124308 neurons 2.9 × 10⁻⁴ → 0.13; GSE142670 neurons 0.37 → not lower; GSE142670 NPCs 2.7 × 10⁻⁵ → 3.3 × 10⁻³. Adjusting for *FOXJ1* alone weakened the signal in every cohort but left it significant in two. The motile-cilia score is built from the same genes as the test set, so adjusting for it removes the set's shared variation by construction. It is shown only to demonstrate that the signal is a coordinated, programme-level shift with no residual gene-specific effect, and not as independent evidence. The primary-cilia set contains intraflagellar-transport and transition-zone genes that every cilium requires (Fig. S4). GSE67528 NPCs were not included in this analysis.

---

## Supplementary Figure S1. Principal component analysis of each cohort and stage.
PCA of log-CPM values (TMM-normalised; 500 most variable genes). Colour, diagnosis; lines join samples (clones or replicates) from the same donor; in (B), triangles mark day 35 and circles day 135. (A) GSE67528 neurons; (B) GSE124308 neurons; (C) GSE142670 neurons; (D) GSE67528 NPCs; (E) GSE142670 NPCs. Diagnosis did not clearly separate on PC1 or PC2 in any cohort, and clones from one donor were often far apart, so clone- or differentiation-level variation was large relative to the diagnosis effect. In (B), PC1 (35.5%) separates the two time points. On PC2, the fibroblast-derived control donors 101 and 102 lie apart from all other samples; every ASD donor was derived from peripheral blood mononuclear cells, as were controls 114 and 574, and the third fibroblast-derived control (1021) did not separate. A per-sample motile-cilia score correlated with a leading component in several cohorts (GSE124308 PC1, r = 0.84; GSE142670 NPC PC2, r = −0.82; GSE67528 neuron PC1, r = 0.52), so ciliated-cell content is a major source of variation between cultures.

## Supplementary Figure S2. Expression changes of the 40 genes with the strongest pooled neuron signal.
Genes with the same direction of change in all three neuron cohorts, ranked by weighted Stouffer z. Fill, log₂ fold change (ASD vs control), capped at ±2; grey, gene not tested in that cohort (removed by expression filtering). Asterisks, nominal p < 0.05 in that cohort; no gene passed FDR < 0.05 in any cohort or in the meta-analysis. Orange squares, cilia-annotated genes (curated motile set or GO cilium organisation/movement/axoneme assembly). NPC columns are shown for comparison and were not used for selection. Because 33 of the 40 genes reached p < 0.05 in GSE124308 (vs 19 in GSE67528 and 5 in GSE142670), the selection is weighted toward GSE124308, and the figure overstates cross-cohort consistency at the level of individual genes.

## Supplementary Figure S3. Overlap of top-ranked genes between neuron cohorts.
UpSet plot of the 500 genes with the smallest p value in each neuron cohort, counted as overlapping only when the direction of change also matched (universe, 14,574 genes tested in all three). Expected counts are from 2,000 permutations that preserved each list's up/down ratio. Overlaps: GSE67528 ∩ GSE124308, 10 observed vs 8.4 expected (permutation p = 0.34); GSE67528 ∩ GSE142670, 8 vs 8.6 (p = 0.62); GSE124308 ∩ GSE142670, 17 vs 9.1 (p = 0.0085); all three, 0 vs 0.15. Top-ranked genes rarely coincide between cohorts, consistent with a diffuse shift spread across many genes (Fig. 2) rather than a few large-effect genes. [The GSE124308 ∩ GSE142670 excess is exploratory.]

## Supplementary Figure S4. STRING interaction network of cilia leading-edge genes (exploratory).
The 60 cilia leading-edge genes with the most negative pooled neuron z (all lower in ASD), queried in STRING (Homo sapiens, combined score ≥ 0.7). 37 genes with at least one interaction are shown (219 edges). Node size, |pooled z|; orange, motile-cilia-specific genes; green, intraflagellar-transport and transition-zone genes required by all cilia, motile and primary; grey, other cilia-annotated genes. The dense cluster corresponds to axonemal dynein, radial-spoke and CFAP proteins. The network illustrates known protein complexes and was not used for inference.

## Supplementary Figure S5. Gene-level random-effects meta-analysis of the three neuron cohorts.
Pooled log₂ fold change (REML random effects, Hartung–Knapp with ad hoc correction) against −log₁₀ p for 14,574 genes. No gene reached FDR < 0.05; the smallest p was [value from the 07b log]. With three cohorts the Hartung–Knapp test has two degrees of freedom, so gene-level inference has very low power. For this reason the primary test of H1 is direction agreement (Fig. 2), not gene-level significance. Orange, curated motile-cilia genes (41 of 42 tested have a pooled estimate below zero); the 12 with the smallest p are labelled. Motile-cilia genes are co-regulated, so this count is not an independent test (see Fig. 3 and 4 for cameraPR).
