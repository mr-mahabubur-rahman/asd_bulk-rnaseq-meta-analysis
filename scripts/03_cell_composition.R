# ==============================================================================
# 03_cell_composition.R — DOES CULTURE COMPOSITION DRIVE THE SHARED ASD SIGNAL?
# ------------------------------------------------------------------------------
# Protocol Section 11, "Cell composition" sensitivity analysis.
# Motivation: mixed neuronal cultures agreed in direction across 3 cohorts,
# but GFP-sorted pure neurons (GSE142670) pointed the opposite way.
#
#   1. Marker-gene scores per sample: neuron, progenitor, proliferation, glia,
#      non-neural (mesenchymal), pluripotency, glutamatergic, GABAergic.
#      Score = mean z-score of expressed markers (logCPM), z computed within
#      each dataset-stage (and within time point for GSE124308).
#   2. ASD vs control at DONOR level (donor means; t-test and Wilcoxon), so
#      clones/time points never count as extra people.
#   3. Re-fit dream with the neuron score as a covariate and repeat the
#      leave-one-cohort-out replication test (neurons) and the NPC comparison,
#      side by side with the unadjusted results.
#
# Input:  data/raw_counts_primary/*, data/metadata/*_metadata_full.csv,
#         data/annotation/ensembl115_genes.tsv, results/tables/DE_dream/*
# Output: results/tables/composition/*.csv
#         results/tables/DE_dream_adjNeuron/*.csv
#         results/figures/composition/composition_scores.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
})

set.seed(2026)
param   <- MulticoreParam(workers = 4)
tab_dir <- "results/tables/composition"
fig_dir <- "results/figures/composition"
adj_dir <- "results/tables/DE_dream_adjNeuron"
for (d in c(tab_dir, fig_dir, adj_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- Marker sets (gene symbols) -----------------------------------------------
markers <- list(
  neuron        = c("MAP2", "STMN2", "SYP", "SNAP25", "RBFOX3", "SYT1", "GAP43",
                    "DCX", "TUBB3", "ELAVL4", "INA", "NEFL", "NEFM"),
  progenitor    = c("NES", "SOX2", "PAX6", "VIM", "HES5", "FABP7", "HES1", "NOTCH1"),
  proliferation = c("MKI67", "TOP2A", "CENPF", "CCNB1", "BIRC5", "PCNA", "MCM2"),
  glia          = c("GFAP", "AQP4", "S100B", "SLC1A3", "ALDH1L1", "OLIG2"),
  non_neural    = c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "TAGLN", "DCN", "LUM", "POSTN"),
  pluripotency  = c("POU5F1", "NANOG", "LIN28A", "DPPA4"),
  glutamatergic = c("SLC17A7", "SLC17A6", "GRIN2B", "TBR1", "NEUROD6"),
  gabaergic     = c("GAD1", "GAD2", "SLC32A1", "DLX5", "DLX6")
)

genes <- read_tsv("data/annotation/ensembl115_genes.tsv",
                  col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE)
sym_of <- setNames(genes$external_gene_name, genes$ensembl_gene_id)

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
           donor = factor(donor))
}

# --- Analyses: dataset, stage filter, base model --------------------------------
analyses <- tribble(
  ~label,              ~gse,        ~stage,   ~base,
  "GSE67528_Neuron",   "GSE67528",  "Neuron", "~ condition + (1|donor)",
  "GSE124308_both",    "GSE124308", NA,       "~ condition + timepoint + (1|donor)",
  "GSE142670_Neuron",  "GSE142670", "Neuron", "~ condition + (1|donor)",
  "GSE67528_NPC",      "GSE67528",  "NPC",    "~ condition + (1|donor)",
  "GSE142670_NPC",     "GSE142670", "NPC",    "~ condition + (1|donor)"
)

# --- 1. Scores ----------------------------------------------------------------
score_one <- function(label, gse, stage) {
  counts <- read_counts(gse)
  meta   <- read_meta(gse)
  if (!is.na(stage)) meta <- filter(meta, stage == !!stage)
  meta   <- filter(meta, sample %in% colnames(counts))
  block  <- if ("timepoint" %in% names(meta)) meta$timepoint else rep("all", nrow(meta))

  dge <- calcNormFactors(DGEList(counts[, meta$sample]))
  lcpm <- cpm(dge, log = TRUE, prior.count = 2)
  sym  <- sym_of[rownames(lcpm)]

  zscore <- function(v) ave(v, block, FUN = function(x) {
    s <- sd(x); if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s })

  sc <- imap_dfc(markers, function(g, nm) {
    idx <- which(sym %in% g & rowMeans(lcpm) > 1)      # expressed markers only
    if (length(idx) < 2) return(tibble(!!nm := NA_real_))
    z <- t(apply(lcpm[idx, , drop = FALSE], 1, zscore))
    tibble(!!nm := colMeans(z))
  })
  bind_cols(meta |> select(sample, donor, condition, any_of("timepoint")), sc) |>
    mutate(analysis = label, .before = 1)
}

scores <- pmap_dfr(analyses |> select(label, gse, stage), score_one)
write_csv(scores, file.path(tab_dir, "composition_scores_by_sample.csv"))

# --- 2. Donor-level ASD vs control ---------------------------------------------
donor_means <- scores |>
  group_by(analysis, donor, condition) |>
  summarise(across(all_of(names(markers)), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

tests <- donor_means |>
  pivot_longer(all_of(names(markers)), names_to = "score", values_to = "value") |>
  filter(is.finite(value)) |>
  group_by(analysis, score) |>
  summarise(
    n_asd = sum(condition == "ASD"), n_ctrl = sum(condition == "Control"),
    diff_asd_minus_ctrl = mean(value[condition == "ASD"]) - mean(value[condition == "Control"]),
    t_p = tryCatch(t.test(value ~ condition)$p.value, error = function(e) NA_real_),
    wilcox_p = tryCatch(wilcox.test(value ~ condition, exact = FALSE)$p.value,
                        error = function(e) NA_real_),
    .groups = "drop") |>
  arrange(analysis, t_p)
message("\n=== Composition: ASD vs control (donor means; diff in SD units) ===")
print(tests |> mutate(across(c(diff_asd_minus_ctrl, t_p, wilcox_p), \(x) signif(x, 2))),
      n = Inf, width = Inf)
write_csv(tests, file.path(tab_dir, "composition_tests_donor_level.csv"))

p <- donor_means |>
  pivot_longer(all_of(names(markers)), names_to = "score", values_to = "value") |>
  filter(score %in% c("neuron", "progenitor", "proliferation", "non_neural")) |>
  ggplot(aes(condition, value, colour = condition)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50") +
  geom_jitter(width = 0.15, size = 1.8) +
  facet_grid(score ~ analysis, scales = "free_y") +
  scale_colour_manual(values = c(Control = "#2166ac", ASD = "#b2182b"), guide = "none") +
  labs(x = NULL, y = "Marker score (mean z, donor mean)",
       title = "Culture composition by diagnosis, one point per donor") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "composition_scores.png"), p, width = 11, height = 8,
       dpi = 300, bg = "white")

# --- 3. Re-fit dream with the neuron score ---------------------------------------
run_dream_adj <- function(label, gse, stage, base) {
  counts <- read_counts(gse)
  meta <- read_meta(gse)
  if (!is.na(stage)) meta <- filter(meta, stage == !!stage)
  meta <- meta |> filter(sample %in% colnames(counts)) |>
    left_join(scores |> filter(analysis == label) |> select(sample, neuron), by = "sample") |>
    as.data.frame()
  if ("timepoint" %in% names(meta)) meta$timepoint <- factor(meta$timepoint)
  rownames(meta) <- meta$sample
  form <- as.formula(paste(base, "+ neuron"))
  message("\n--- ", label, " | ", deparse(form))

  dge  <- DGEList(counts[, meta$sample])
  keep <- filterByExpr(dge, group = meta$condition)
  dge  <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  vobj <- voomWithDreamWeights(dge, form, meta, BPPARAM = param)
  fit  <- eBayes(dream(vobj, form, meta, BPPARAM = param))
  tt   <- topTable(fit, coef = "conditionASD", number = Inf, sort.by = "none")
  res  <- tibble(Gene_ID = rownames(tt), log2FoldChange = tt$logFC,
                 lfcSE = abs(tt$logFC / tt$t), stat = tt$t,
                 pvalue = tt$P.Value, padj = tt$adj.P.Val,
                 Gene_Symbol = unname(sym_of[rownames(tt)]))
  write_csv(res, file.path(adj_dir, paste0(label, ".csv")))
  message("Genes at FDR < 0.05: ", sum(res$padj < 0.05, na.rm = TRUE))
  invisible(res)
}
pwalk(analyses, run_dream_adj)

# --- 4. Replication: unadjusted vs neuron-adjusted --------------------------------
donors <- read_csv("results/tables/DE_dream/summary.csv", show_col_types = FALSE) |>
  select(analysis, donors) |> deframe()
stouffer_z <- function(lfc, p, w) {
  z <- sign(lfc) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
  sum(w * z) / sqrt(sum(w^2))
}
load_de <- function(dir, a) read_csv(file.path(dir, paste0(a, ".csv")),
                                     show_col_types = FALSE) |>
  select(Gene_ID, log2FoldChange, pvalue) |> mutate(dataset = a)

replicate_loo <- function(dir, sets, version) {
  long <- map_dfr(sets, ~ load_de(dir, .x))
  map_dfr(sets, function(held) {
    top <- long |> filter(dataset != held) |>
      group_by(Gene_ID) |> filter(n() == length(sets) - 1) |>
      summarise(z = stouffer_z(log2FoldChange, pvalue, sqrt(donors[dataset])),
                .groups = "drop") |>
      arrange(desc(abs(z))) |> head(500)
    h <- long |> filter(dataset == held) |> inner_join(top, by = "Gene_ID")
    k <- sum(sign(h$log2FoldChange) == sign(h$z))
    tibble(version = version, held_out = held, genes = nrow(h), share = k / nrow(h),
           binom_p_two_sided = binom.test(k, nrow(h), 0.5)$p.value)
  })
}
neuron_sets <- c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron")
rep_tbl <- bind_rows(
  replicate_loo("results/tables/DE_dream", neuron_sets, "unadjusted"),
  replicate_loo(adj_dir, neuron_sets, "neuron-adjusted"))
message("\n=== Neurons: leave-one-cohort-out direction agreement (top 500) ===")
print(rep_tbl, width = Inf)
write_csv(rep_tbl, file.path(tab_dir, "replication_unadj_vs_adjNeuron.csv"))

npc_tbl <- bind_rows(
  replicate_loo("results/tables/DE_dream", c("GSE67528_NPC", "GSE142670_NPC"), "unadjusted"),
  replicate_loo(adj_dir, c("GSE67528_NPC", "GSE142670_NPC"), "neuron-adjusted"))
message("\n=== NPC: direction agreement between the two cohorts (top 500) ===")
print(npc_tbl, width = Inf)
write_csv(npc_tbl, file.path(tab_dir, "npc_agreement_unadj_vs_adjNeuron.csv"))
