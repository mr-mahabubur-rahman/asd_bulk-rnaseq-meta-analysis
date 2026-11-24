# ==============================================================================
# 02_meta_analysis.R — STAGE-MATCHED META-ANALYSIS (protocol Sections 10-11)
# ------------------------------------------------------------------------------
# Neurons (primary): GSE67528_Neuron, GSE124308_both, GSE142670_Neuron
#   - per-gene random-effects model, REML tau^2, Hartung-Knapp with ad hoc
#     correction (metafor, test = "adhoc"; protocol amendment 1)
#   - genes tested in >= 3 datasets; BH FDR across genes
#   - robust gene = FDR < 0.05 AND pooled p < 0.05 with any one dataset left out
#   - complementary weighted Stouffer (weights = sqrt(donors))
#   - leave-one-cohort-out replication: direction agreement + pi1 in held-out set
#   - sensitivity: GSE142670 sorted neurons instead of mixed neurons
# NPC (secondary): only 2 datasets (GSE67528_NPC, GSE142670_NPC) -> no
#   meta-analysis (protocol: < 3 datasets); direct comparison only.
#
# NOTE: GSE142670 currently uses the AUTHORS' counts (GENCODE 19, HTSeq);
# final results use uniform Salmon re-quantification.
#
# Input:  results/tables/DE_dream/<analysis>.csv, summary.csv (donor counts)
# Output: results/tables/meta/*.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(limma)   # propTrueNull()
})

in_dir  <- "results/tables/DE_dream"
out_dir <- "results/tables/meta"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

donors <- read_csv(file.path(in_dir, "summary.csv"), show_col_types = FALSE) |>
  select(analysis, donors) |> deframe()

read_de <- function(a) {
  read_csv(file.path(in_dir, paste0(a, ".csv")), show_col_types = FALSE) |>
    select(Gene_ID, Gene_Symbol, log2FoldChange, lfcSE, stat, pvalue) |>
    filter(is.finite(log2FoldChange), is.finite(lfcSE), lfcSE > 0) |>
    mutate(dataset = a)
}

# --- Per-gene random-effects model --------------------------------------------
# Hartung-Knapp with the ad hoc variance correction (Roever et al. 2015):
# the HK variance is never allowed below the standard random-effects variance.
# Plain "knha" gives artificially narrow CIs when, with k = 3, the estimates
# happen to agree closely (protocol amendment 1).
fit_rem <- function(yi, sei) {
  fit <- tryCatch(
    suppressWarnings(rma(yi = yi, sei = sei, method = "REML", test = "adhoc")),
    error = function(e) tryCatch(
      suppressWarnings(rma(yi = yi, sei = sei, method = "DL", test = "adhoc")),
      error = function(e) NULL))
  if (is.null(fit)) return(c(est = NA, se = NA, ci_lb = NA, ci_ub = NA, p = NA,
                             tau2 = NA, I2 = NA, Q = NA))
  c(est = fit$b[1], se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
    p = fit$pval, tau2 = fit$tau2, I2 = fit$I2, Q = fit$QE)
}

stouffer <- function(lfc, p, w) {
  z <- sign(lfc) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
  zc <- sum(w * z) / sqrt(sum(w^2))
  c(z = zc, p = 2 * pnorm(-abs(zc)))
}

run_meta <- function(sets, label, min_k = 3) {
  message("\n==================== ", label, " ====================")
  message("Datasets: ", paste(sets, collapse = ", "))
  long <- map_dfr(sets, read_de)
  k_tbl <- long |> count(Gene_ID, name = "k")
  genes <- k_tbl |> filter(k >= min_k) |> pull(Gene_ID)
  message("Genes in >= ", min_k, " datasets: ", length(genes))
  long <- long |> filter(Gene_ID %in% genes)

  # Random effects, gene by gene
  message("Fitting random-effects models (a few minutes)...")
  rem <- long |>
    group_by(Gene_ID) |>
    summarise(Gene_Symbol = first(na.omit(Gene_Symbol)),
              k = n(),
              same_sign = abs(sum(sign(log2FoldChange))) == n(),
              rem = list(fit_rem(log2FoldChange, lfcSE)),
              st  = list(stouffer(log2FoldChange, pvalue, sqrt(donors[dataset]))),
              .groups = "drop") |>
    unnest_wider(rem) |>
    unnest_wider(st, names_sep = "_") |>
    mutate(fdr = p.adjust(p, "BH"),
           stouffer_fdr = p.adjust(st_p, "BH"))

  # Leave-one-dataset-out for genes near the threshold (FDR < 0.10)
  cand <- rem |> filter(fdr < 0.10) |> pull(Gene_ID)
  loo <- map_dfr(sets, function(drop) {
    long |> filter(Gene_ID %in% cand, dataset != drop) |>
      group_by(Gene_ID) |>
      summarise(p_loo = fit_rem(log2FoldChange, lfcSE)[["p"]], .groups = "drop") |>
      mutate(dropped = drop)
  })
  loo_max <- loo |> group_by(Gene_ID) |>
    summarise(max_p_leave_one_out = suppressWarnings(max(p_loo, na.rm = TRUE)),
              .groups = "drop")
  rem <- rem |> left_join(loo_max, by = "Gene_ID") |>
    mutate(robust = fdr < 0.05 & !is.na(max_p_leave_one_out) &
             max_p_leave_one_out < 0.05) |>
    arrange(p)

  write_csv(rem, file.path(out_dir, paste0(label, "_REM.csv")))

  summ <- tibble(
    analysis = label, datasets = length(sets), genes_tested = nrow(rem),
    pi1 = 1 - propTrueNull(rem$p[!is.na(rem$p)], method = "convest"),
    fdr_lt_0.05 = sum(rem$fdr < 0.05, na.rm = TRUE),
    robust = sum(rem$robust, na.rm = TRUE),
    stouffer_fdr_lt_0.05 = sum(rem$stouffer_fdr < 0.05, na.rm = TRUE),
    min_fdr = min(rem$fdr, na.rm = TRUE))
  message("\nSummary:")
  print(summ, width = Inf)
  message("\nTop 20 genes by random-effects p:")
  print(rem |> select(Gene_Symbol, k, est, ci_lb, ci_ub, p, fdr, I2,
                      same_sign, st_p, robust) |> head(20), width = Inf)

  # Leave-one-cohort-out replication
  message("\nLeave-one-cohort-out replication:")
  rep_tbl <- map_dfr(sets, function(held) {
    others <- long |> filter(dataset != held) |>
      group_by(Gene_ID) |>
      summarise(st = list(stouffer(log2FoldChange, pvalue, sqrt(donors[dataset]))),
                lfc = mean(log2FoldChange), .groups = "drop") |>
      unnest_wider(st, names_sep = "_") |>
      mutate(fdr = p.adjust(st_p, "BH"))
    called <- others |> filter(fdr < 0.05)
    rule <- "FDR < 0.05 in the other datasets"
    if (nrow(called) < 50) {
      called <- others |> arrange(st_p) |> head(500)
      rule <- "top 500 by p in the other datasets (fewer than 50 at FDR < 0.05)"
    }
    h <- long |> filter(dataset == held) |>
      inner_join(called |> select(Gene_ID, z_other = st_z), by = "Gene_ID")
    k <- sum(sign(h$log2FoldChange) == sign(h$z_other))
    tibble(held_out = held, genes = nrow(h), rule = rule,
           same_direction = k, share = k / nrow(h),
           binom_p = binom.test(k, nrow(h), 0.5, alternative = "greater")$p.value,
           pi1_held_out = 1 - propTrueNull(h$pvalue, method = "convest"))
  })
  print(rep_tbl, width = Inf)
  write_csv(rep_tbl, file.path(out_dir, paste0(label, "_replication.csv")))

  list(summary = summ, replication = rep_tbl)
}

# --- Neurons: primary + sensitivity --------------------------------------------
neuron <- run_meta(c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron"),
                   "neuron_primary")
neuron_sorted <- run_meta(c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron_sorted"),
                          "neuron_sens_sorted")

# --- NPC: direct comparison only (2 datasets) ----------------------------------
message("\n==================== NPC: direct comparison ====================")
a <- read_de("GSE67528_NPC"); b <- read_de("GSE142670_NPC")
m <- inner_join(a, b, by = "Gene_ID", suffix = c("_a", "_b"))
npc_tbl <- map_dfr(c(100, 250, 500), function(n) {
  bind_rows(
    { t <- m |> arrange(pvalue_a) |> head(n)
      k <- sum(sign(t$log2FoldChange_a) == sign(t$log2FoldChange_b))
      tibble(discovery = "GSE67528_NPC", top_n = n, share = k / n,
             binom_p = binom.test(k, n, 0.5, alternative = "greater")$p.value) },
    { t <- m |> arrange(pvalue_b) |> head(n)
      k <- sum(sign(t$log2FoldChange_a) == sign(t$log2FoldChange_b))
      tibble(discovery = "GSE142670_NPC", top_n = n, share = k / n,
             binom_p = binom.test(k, n, 0.5, alternative = "greater")$p.value) })
}) |> mutate(spearman_t = cor(m$stat_a, m$stat_b, method = "spearman"))
print(npc_tbl, width = Inf)
write_csv(npc_tbl, file.path(out_dir, "npc_direct_comparison.csv"))

write_csv(bind_rows(neuron$summary, neuron_sorted$summary),
          file.path(out_dir, "meta_summary.csv"))
