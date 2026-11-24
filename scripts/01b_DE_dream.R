# ==============================================================================
# 01b_DE_dream.R — DONOR-AWARE DIFFERENTIAL EXPRESSION (replaces 01_DESeq2.R)
# ------------------------------------------------------------------------------
# Why: samples are not independent people. GSE67528 = 14 donors x 3 stages
# (1-4 clones per donor per stage); GSE124308 = 11 donors x 2 time points.
# dream (variancePartition) fits a linear mixed model with a random effect
# for donor, so repeated samples from one person are not counted as extra
# individuals. Counts: primary-assembly genes only (from 01a).
#
# Analyses
#   GSE67528  : one model per stage (iPSC, NPC, Neuron)
#               ~ condition + (1|donor)
#   GSE142670 : one model per stage (NPC, Neuron, Neuron_sorted)
#               ~ condition + (1|donor)      [authors' counts for now]
#   GSE124308 : main       ~ condition + timepoint + (1|donor)
#               per time   ~ condition            (d35, d135; 1 sample/donor)
#               sensitivity~ condition + timepoint + ipsc_source + (1|donor)
#               (all 6 ASD lines are PBMC/Sendai; 3 of 5 controls are
#                fibroblast-derived, so source is partly confounded)
#
# Output: results/tables/DE_dream/<analysis>.csv
#         columns match the old DESeq2 output (Gene_ID, log2FoldChange,
#         lfcSE, stat, pvalue, padj) so the meta-analysis script still works
#         results/tables/DE_dream/summary.csv
# ==============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)

set.seed(2026)
param   <- MulticoreParam(workers = 4)   # lower to 2 if the laptop struggles
out_dir <- "results/tables/DE_dream"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Inputs -------------------------------------------------------------------
genes <- read_tsv("data/annotation/ensembl115_genes.tsv",
                  col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE) |>
  select(Gene_ID = ensembl_gene_id, Gene_Symbol = external_gene_name,
         Gene_Biotype = gene_biotype)

read_counts <- function(geo_id) {
  x <- read_csv(paste0("data/raw_counts_primary/", geo_id, "_raw_counts.csv"),
                show_col_types = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  m
}

read_meta <- function(geo_id) {
  read_csv(paste0("data/metadata/", geo_id, "_metadata_full.csv"),
           show_col_types = FALSE) |>
    mutate(condition = factor(condition, levels = c("Control", "ASD")),
           donor     = factor(donor)) |>
    as.data.frame()
}

# --- Core function ------------------------------------------------------------
# counts: genes x samples matrix; meta: data.frame with one row per sample
run_dream <- function(counts, meta, formula, label, coef = "conditionASD") {
  message("\n--- ", label, " ---")
  no_counts <- setdiff(meta$sample, colnames(counts))
  if (length(no_counts) > 0) {
    message("Samples without counts (excluded): ", paste(no_counts, collapse = ", "))
    meta <- meta[meta$sample %in% colnames(counts), ]
  }
  meta   <- droplevels(meta)
  rownames(meta) <- meta$sample
  counts <- counts[, meta$sample, drop = FALSE]
  message("Samples: ", ncol(counts),
          " | donors: ", n_distinct(meta$donor),
          " (", paste(names(table(distinct(meta, donor, condition)$condition)),
                      table(distinct(meta, donor, condition)$condition),
                      collapse = ", "), ")")

  # Filter + TMM normalisation (works with non-integer RSEM counts too)
  dge  <- DGEList(counts)
  keep <- filterByExpr(dge, group = meta$condition)
  dge  <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  message("Genes after filterByExpr: ", nrow(dge))

  # voom weights + mixed model
  vobj <- voomWithDreamWeights(dge, formula, meta, BPPARAM = param)
  fit  <- dream(vobj, formula, meta, BPPARAM = param)
  fit  <- eBayes(fit)

  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  res <- tibble(
    Gene_ID        = rownames(tt),
    log2FoldChange = tt$logFC,
    lfcSE          = abs(tt$logFC / tt$t),   # SE implied by moderated t
    stat           = tt$t,
    pvalue         = tt$P.Value,
    padj           = tt$adj.P.Val,
    AveExpr        = tt$AveExpr
  ) |>
    left_join(genes, by = "Gene_ID")

  write_csv(res, file.path(out_dir, paste0(label, ".csv")))

  tibble(analysis = label,
         samples  = ncol(counts),
         donors   = n_distinct(meta$donor),
         genes_tested = nrow(res),
         padj_lt_0.05 = sum(res$padj < 0.05, na.rm = TRUE),
         padj_lt_0.05_lfc_gt_1 = sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1,
                                     na.rm = TRUE),
         min_padj = min(res$padj, na.rm = TRUE))
}

# Datasets to run: from the command line, e.g.
#   pixi run Rscript scripts/01b_DE_dream.R GSE142670
# default = all datasets with a block below
run_ids <- commandArgs(trailingOnly = TRUE)
if (length(run_ids) == 0) run_ids <- c("GSE67528", "GSE124308", "GSE142670")

summaries <- list()

# --- GSE67528: per stage, donor random effect ---------------------------------
if ("GSE67528" %in% run_ids) {
c1 <- read_counts("GSE67528")
m1 <- read_meta("GSE67528") |>
  mutate(stage = factor(stage, levels = c("iPSC", "NPC", "Neuron")))

for (st in levels(m1$stage)) {
  summaries[[paste0("GSE67528_", st)]] <- run_dream(
    c1, filter(m1, stage == st),
    formula = ~ condition + (1 | donor),
    label   = paste0("GSE67528_", st)
  )
}

}

# --- GSE124308: neurons at d35 + d135 -----------------------------------------
if ("GSE124308" %in% run_ids) {
c2 <- read_counts("GSE124308")
m2 <- read_meta("GSE124308") |>
  mutate(timepoint   = factor(timepoint, levels = c("d35", "d135")),
         ipsc_source = factor(ipsc_source))

# Main: ASD effect averaged over time, donor random effect
summaries[["GSE124308_both"]] <- run_dream(
  c2, m2, formula = ~ condition + timepoint + (1 | donor),
  label = "GSE124308_both")

# Per time point (one sample per donor, so no random effect needed)
for (tp in levels(m2$timepoint)) {
  summaries[[paste0("GSE124308_", tp)]] <- run_dream(
    c2, filter(m2, timepoint == tp), formula = ~ condition,
    label = paste0("GSE124308_", tp))
}

# Sensitivity: adjust for iPSC source (fibroblast vs PBMC)
summaries[["GSE124308_both_adjSource"]] <- run_dream(
  c2, m2, formula = ~ condition + timepoint + ipsc_source + (1 | donor),
  label = "GSE124308_both_adjSource")

}

# --- GSE142670: NPC, mixed neurons, sorted neurons (authors' counts) -----------
# 12 donors (6 ASD / 6 control), 1-3 clones per donor per stage
if ("GSE142670" %in% run_ids) {
c3 <- read_counts("GSE142670")
m3 <- read_meta("GSE142670") |>
  mutate(stage = factor(stage, levels = c("NPC", "Neuron", "Neuron_sorted")))

for (st in levels(m3$stage)) {
  summaries[[paste0("GSE142670_", st)]] <- run_dream(
    c3, filter(m3, stage == st),
    formula = ~ condition + (1 | donor),
    label   = paste0("GSE142670_", st)
  )
}
}

# --- Summary (rows for re-run analyses are replaced, others kept) --------------
summary_tbl <- bind_rows(summaries)
summary_file <- file.path(out_dir, "summary.csv")
if (file.exists(summary_file)) {
  old_tbl <- read_csv(summary_file, show_col_types = FALSE) |>
    filter(!analysis %in% summary_tbl$analysis)
  summary_tbl <- bind_rows(old_tbl, summary_tbl) |> arrange(analysis)
}
write_csv(summary_tbl, summary_file)
message("\n==================== SUMMARY ====================")
print(summary_tbl, width = Inf)
