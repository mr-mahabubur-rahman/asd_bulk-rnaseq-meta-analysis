# Transcriptomic alterations in iPSC-derived neural models of idiopathic autism

**A systematic review and donor-aware meta-analysis of re-processed RNA-seq data: protocol, analysis pipeline and discovery-stage results**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Status (24 Nov 2026)**
> - **Protocol v1.0 is frozen** ([`docs/protocol_v1.0.pdf`](docs/protocol_v1.0.pdf)). It was not prospectively registered. Its git commits serve as timestamps: v0.3 is commit `3c2e0a9` (24 Nov 2026, 13:24 UTC+6) and v1.0 is commit `33d4e50` (17:03 UTC+6).
> - **The discovery analysis is complete** for three cohorts (GSE67528, GSE124308, GSE142670). Its results generated the hypotheses H1–H4. They are **not** confirmatory.
> - **Confirmation is pending.** It requires a systematic search, uniform Salmon re-quantification of all eligible datasets, and tests of H1–H4 in datasets not used for discovery.
>
> This repository **replaces** an earlier pilot analysis that reported 39 meta-analysis DEGs and an interferon/MHC-I signature. Those results are **withdrawn**: see [What changed from the pilot](#what-changed-from-the-pilot).

---

## Authors

| | Name | Affiliation | ORCID |
|---|---|---|---|
| 1 | **Mahabubur Rahman** | Molecular Diagnostics and NGS Laboratory, BIRDEM General Hospital, Dhaka, Bangladesh | [0000-0002-5535-5816](https://orcid.org/0000-0002-5535-5816) |
| 2 | **Mofizul Islam** | Department of Neurology, BIRDEM General Hospital, Dhaka, Bangladesh | — |

Correspondence: Mahabubur Rahman [email to be added].

---

## Question

Do neural progenitor cells (NPCs) and neurons derived from induced pluripotent stem cells (iPSCs) of people with **idiopathic** autism show reproducible transcriptomic differences from controls once **donor** (not sample) is treated as the unit of replication?

## Pre-specified hypotheses (protocol v1.0)

| | Hypothesis | Primary test (confirmation set) |
|---|---|---|
| **H1** | ASD-associated expression changes replicate across cohorts in direction | Binomial test of direction agreement for the discovery top-500 genes |
| **H2** | A FOXJ1-associated ciliogenesis programme is lower in ASD, reflecting fewer or less mature FOXJ1-positive ciliated cells rather than loss of primary cilia from neurons | cameraPR, curated motile-cilia set, one-sided (lower) |
| **H3** | Hallmark MYC targets V1 are higher in ASD NPCs | cameraPR, one-sided (higher) |
| **H4** | Hallmark interferon-α response is lower in ASD | cameraPR, one-sided (lower) |

The four primary tests are Holm-corrected. Results from multiple confirmation datasets are combined by weighted Stouffer z (weights = √donors).

---

## Datasets

| Dataset | Reference | Model and stages | ASD / control donors | Role |
|---|---|---|---|---|
| GSE67528 | Marchetto et al., *Mol Psychiatry* 2017 | 2D: iPSC, NPC, neurons | 8 / 6 | Discovery |
| GSE124308 | DeRosa et al., *Sci Rep* 2018 | 2D cortical neurons, day 35 and day 135 | 6 / 5 | Discovery (iPSC source partly confounded with diagnosis) |
| GSE142670 | Griesi-Oliveira et al., *Mol Psychiatry* 2021 | 2D: NPC, mixed neurons, FACS-sorted neurons | 6 / 6 | Discovery |
| GSE61476 | Mariani et al., *Cell* 2015 | Organoids, day 0, 11, 31 | 4 / 8 (family design) | Confirmation (organoid analysis) |
| syn8118403 | Adhya et al., *Biol Psychiatry* 2021 | 2D cortical neurons, day 35 | 3 / 3 | Confirmation, if accessible |
| E-MTAB-6018 | Schafer et al., *Nat Neurosci* 2019 | 2D neural time course | 8 / 5 | Donor-overlap check against GSE67528 first |

Further candidates and exclusions are listed in the protocol (Study records). The final dataset list will come from the systematic search.

---

## Discovery-stage findings (hypothesis-generating)

| Finding | Result |
|---|---|
| Single-gene differential expression | **No gene at FDR < 0.05** in any cohort, or in the 3-cohort random-effects meta-analysis (14,574 genes; smallest p = 0.045) |
| **Direction replication (neurons)** | Leave-one-cohort-out agreement for the top 500 genes: 61–63% (one-sided binomial p ≤ 9 × 10⁻⁷); 58–69% after neuron-score adjustment |
| Direction replication (NPCs) | 56–63% (p = 0.02 to 1 × 10⁻⁷) |
| FACS-sorted GSE142670 neurons | 33% agreement, i.e. **reversed** (two-sided p = 5 × 10⁻¹⁵), from only 3 ASD / 3 control donors. Most replicating genes that are lower in ASD are removed by sorting (non-neuronal), so the reversal most likely reflects sorting purity, not opposite neuronal biology |
| Cell composition | No difference by diagnosis in neuronal or progenitor marker scores |
| **Cilia (motile programme)** | Lower in ASD in GSE67528 neurons, GSE124308 neurons and GSE142670 NPCs (cameraPR); not in GSE142670 neurons. Adjusting for *FOXJ1* weakens but does not remove the signal |
| **MYC targets / cell cycle (NPCs)** | Higher in ASD NPCs in both cohorts (cameraPR FDR 4.7 × 10⁻³ and 3.9 × 10⁻⁷) |
| Interferon | Lower in ASD in GSE124308 only. Within that cohort it is robust to adjusting for iPSC source and to excluding the two outlying control donors, but it does not appear in GSE67528 or GSE142670 (cohort-specific) |
| Sensitivity: GSE124308 without control donors 101/102 | Direction replication still holds (55–58%, p ≤ 0.009) |

**Interpretation.** The signal is a weak but reproducible shift spread across many genes, not a set of large-effect genes. Much of the part that is lower in ASD sits in genes expressed by non-neuronal cells (they are removed by FACS sorting), and the cilia result most likely reflects how many ciliated (FOXJ1-positive) cells each culture contains. The confirmation stage tests these points in independent data.

Figures: `results/figures/submission/` (Fig 2–4, S1–S5). Draft legends: [`docs/figure_legends_draft.md`](docs/figure_legends_draft.md).

---

## What changed from the pilot

The pilot (archived in `legacy/` and at git tag `pilot-v0`) analysed GSE67528 and GSE124308 with DESeq2 + SVA and reported 39 DEGs, KLF4, and interferon/MHC-I hub genes. Re-analysis showed that those results were artefacts:

| Problem in the pilot | Fix in this version |
|---|---|
| **Alternative-haplotype and patch gene IDs** were counted as separate genes. The HLA/MHC "hits" came from genes duplicated across ALT contigs | Restricted to primary-assembly genes (chromosomes 1–22, X, Y, MT; Ensembl release 115) |
| **Pseudoreplication**: clones and replicates from the same donor were treated as independent samples | `dream` mixed models with donor as a random effect |
| **Developmental stages were pooled** (iPSC, NPC and neurons in one model) | Separate analyses for each stage |
| Gene overlap (Venn diagram) plus random-effects meta-analysis of two cohorts | Pre-specified protocol, 3 discovery cohorts, direction-replication test, correlation-aware gene-set tests (cameraPR) |
| Confounding of iPSC source with diagnosis in GSE124308 not addressed | Documented; sensitivity analyses that adjust for source and that exclude the two outlying fibroblast-derived control donors (`09_sensitivity_checks.R`) |

---

## Pipeline

All analysis runs in a [pixi](https://pixi.sh) environment (`pixi.toml`; conda-forge and bioconda, linux-64; developed under WSL Ubuntu).

| Step | Script | What it does |
|---|---|---|
| 00a | `00a_get_ensembl_annotation.sh` | Downloads the Ensembl 115 GTF and extracts gene ID, chromosome, symbol and biotype |
| 00b | `00b_geo_metadata.R` | Downloads full GEO sample annotation (pilot cohorts) |
| 00c | `00c_tidy_metadata.R` | Tidy metadata for GSE67528 and GSE124308: donor, clone, stage, time point, iPSC source |
| 00d | `00d_inspect_geo.R` | Prints GEO fields, processing notes and supplementary files for new accessions |
| 00e | `00e_tidy_metadata_new.R` | Tidy metadata for GSE142670 and GSE61476 |
| 00f | `00f_get_author_counts.R` | Downloads per-sample author count files and assembles a count matrix |
| 01a | `01a_filter_primary_genes.R` | Keeps primary-assembly genes; removes ALT/patch and retired IDs |
| 01b | `01b_DE_dream.R` | Donor-aware differential expression (filterByExpr, TMM, voomWithDreamWeights, `dream`) |
| 01c | `01c_DE_diagnostics.R` | π1, p-value histograms, cross-cohort concordance |
| 02 | `02_meta_analysis.R` | Random-effects meta-analysis (REML, Hartung–Knapp ad hoc), weighted Stouffer, leave-one-cohort-out replication |
| 03 | `03_cell_composition.R` | Marker-gene composition scores, donor-level tests, neuron-adjusted refits |
| 04 | `04_pathway_analysis.R` | fgsea screening (Hallmark, GO:BP, Reactome) |
| 05 | `05_cilia_check.R` | cameraPR, motile vs primary cilia sets, leading edge, ependymal markers |
| 06 | `06_cilia_foxj1_adjusted.R` | Refits adjusted for FOXJ1 or a motile-cilia score |
| 07 | `07_figures.R` | First-draft descriptive figures (superseded by 07b) |
| 07b | `07b_figures_revised.R` | Submission figures: Fig 3, S1–S5 |
| 08 | `08_main_figures.R` | Submission figures: Fig 2 (H1) and Fig 4 (H2) |
| 09 | `09_sensitivity_checks.R` | GSE124308 without donors 101/102; sorted vs mixed GSE142670 neurons |

### How to run

```bash
git clone https://github.com/mr-mahabubur-rahman/asd_bulk-rnaseq-meta-analysis.git
cd asd_bulk-rnaseq-meta-analysis
pixi install

bash scripts/00a_get_ensembl_annotation.sh
pixi run geo-metadata                                   # 00b
pixi run tidy-metadata                                  # 00c
pixi run Rscript scripts/00d_inspect_geo.R GSE142670 GSE61476
pixi run Rscript scripts/00e_tidy_metadata_new.R
pixi run Rscript scripts/00f_get_author_counts.R GSE142670 GSE61476
pixi run Rscript scripts/01a_filter_primary_genes.R GSE67528 GSE124308 GSE142670
pixi run de-dream                                       # 01b (slowest step)
pixi run de-diagnostics                                 # 01c
for s in 02_meta_analysis 03_cell_composition 04_pathway_analysis 05_cilia_check \
         06_cilia_foxj1_adjusted 08_main_figures 07b_figures_revised 09_sensitivity_checks; do
  pixi run bash -c "mkdir -p logs && Rscript scripts/$s.R 2>&1 | tee logs/$s.log"
done
```

Every step writes a log to `logs/`. Scripts 07 and 07b query the STRING API, and 00a/00b/00f need internet access.

### Input data

The discovery analysis uses **deposited count matrices**: the series-level tables for GSE67528 and GSE124308, and the authors' per-sample count files for GSE142670. All three are harmonised to Ensembl 115 primary-assembly gene IDs. The protocol's main analysis will use **uniform Salmon re-quantification from FASTQ** for every dataset. The deposited counts will be kept as a sensitivity analysis.

---

## Repository structure

```
asd_bulk-rnaseq-meta-analysis/
├── README.md
├── LICENSE
├── pixi.toml / pixi.lock
├── docs/
│   ├── protocol_v0.3.pdf, protocol_v1.0.pdf   # protocol versions (frozen)
│   └── figure_legends_draft.md
├── scripts/                                   # 00a–09, run in order (see table)
├── legacy/                                    # withdrawn pilot pipeline (read-only)
├── data/
│   ├── annotation/            # Ensembl 115 gene table
│   ├── metadata/              # GEO annotation + tidy metadata per dataset
│   ├── raw_counts/            # deposited count matrices (pilot cohorts)
│   ├── raw_counts_author/     # assembled author counts (GSE142670)
│   └── raw_counts_primary/    # primary-assembly genes only (analysis input)
├── logs/                      # one log per script run
└── results/
    ├── tables/                # DE_dream/, meta/, composition/, pathways, cilia/, figures/
    └── figures/submission/    # Fig 2–4 and S1–S5 (PDF + PNG)
```

---

## Methods references

- **Differential expression:** edgeR (Robinson et al., 2010); limma/voom (Ritchie et al., 2015; Law et al., 2014); dream and variancePartition (Hoffman & Roussos, 2021; Hoffman & Schadt, 2016)
- **Meta-analysis:** metafor (Viechtbauer, 2010); Hartung–Knapp ad hoc correction (Röver et al., 2015); weighted Z (Whitlock, 2005)
- **Gene sets:** fgsea (Korotkevich et al.); camera (Wu & Smyth, 2012); MSigDB Hallmark (Liberzon et al., 2015)
- **Reporting:** PRISMA-P 2015 (Moher et al.) and PRISMA 2020 (Page et al.)

Full citations are in the protocol.

## Citation

Until the protocol and paper are published, please cite this repository:

> Rahman M, Islam M. Transcriptomic alterations in iPSC-derived neural models of idiopathic autism: protocol for a systematic review and donor-aware meta-analysis of re-processed RNA-seq data. Version 1.0. 2026. https://github.com/mr-mahabubur-rahman/asd_bulk-rnaseq-meta-analysis

## License

MIT — see [LICENSE](LICENSE). The datasets belong to their original authors and are subject to the terms of GEO, ArrayExpress, Synapse and NDA.
