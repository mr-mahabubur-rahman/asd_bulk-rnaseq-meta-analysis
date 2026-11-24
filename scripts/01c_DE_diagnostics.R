# ==============================================================================
# 01c_DE_diagnostics.R — IS THERE ANY ASD SIGNAL BELOW THE FDR THRESHOLD?
# ------------------------------------------------------------------------------
# No gene reached FDR < 0.05 in the donor-aware models (except the
# source-adjusted GSE124308 model). This script asks three questions:
#   1. Do p-values pile up near 0 (pi1 = estimated share of non-null genes)?
#   2. Do the two datasets agree in direction more than chance (50%)?
#   3. How do the 11 source-adjusted "hits" behave in the other analyses?
#
# Input:  results/tables/DE_dream/*.csv   (from 01b_DE_dream.R)
# Output: results/tables/DE_dream/diagnostics_*.csv
#         results/figures/DE_dream/pvalue_histograms.png
#         results/figures/DE_dream/neuron_concordance.png
# ==============================================================================

library(tidyverse)
library(limma)   # propTrueNull()

in_dir  <- "results/tables/DE_dream"
fig_dir <- "results/figures/DE_dream"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(in_dir, pattern = "^GSE.*\\.csv$", full.names = TRUE)
res   <- set_names(map(files, read_csv, show_col_types = FALSE),
                   tools::file_path_sans_ext(basename(files)))

# --- 1. p-value distributions and pi1 -----------------------------------------
pi_tbl <- imap_dfr(res, function(df, nm) {
  p <- df$pvalue[!is.na(df$pvalue)]
  tibble(analysis      = nm,
         genes         = length(p),
         pi1           = 1 - propTrueNull(p, method = "convest"),
         p_lt_0.01     = sum(p < 0.01),
         expected_0.01 = round(0.01 * length(p)))
})
message("\n=== 1. Share of genes with a real ASD effect (pi1) ===")
print(pi_tbl, width = Inf)
write_csv(pi_tbl, file.path(in_dir, "diagnostics_pi1.csv"))

p_hist <- bind_rows(res, .id = "analysis") |>
  ggplot(aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.02), fill = "grey60", colour = "white") +
  facet_wrap(~ analysis, scales = "free_y") +
  labs(x = "p-value (ASD vs control)", y = "Genes",
       title = "A real signal shows as a peak near p = 0") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "pvalue_histograms.png"), p_hist,
       width = 10, height = 7, dpi = 300, bg = "white")

# --- 2. Cross-dataset direction agreement --------------------------------------
concordance <- function(a, b, label_a, label_b, top_n = c(100, 250, 500)) {
  m <- inner_join(a, b, by = "Gene_ID", suffix = c("_a", "_b"))
  rho <- cor(m$stat_a, m$stat_b, method = "spearman", use = "complete.obs")
  one_way <- function(disc, rep_col, from, to) {
    map_dfr(top_n, function(n) {
      top <- m |> arrange(.data[[paste0("pvalue_", disc)]]) |> head(n)
      k   <- sum(sign(top[[paste0("log2FoldChange_", disc)]]) ==
                 sign(top[[paste0("log2FoldChange_", rep_col)]]), na.rm = TRUE)
      tibble(discovery = from, replication = to, top_n = n,
             same_direction = k, share = k / n,
             binom_p = binom.test(k, n, 0.5, alternative = "greater")$p.value)
    })
  }
  out <- bind_rows(one_way("a", "b", label_a, label_b),
                   one_way("b", "a", label_b, label_a)) |>
    mutate(shared_genes = nrow(m), spearman_t = rho)
  list(table = out, merged = m)
}

pairs <- list(
  c("GSE67528_Neuron", "GSE124308_both"),
  c("GSE67528_NPC",    "GSE124308_d35")
)
conc <- map(pairs, ~ concordance(res[[.x[1]]], res[[.x[2]]], .x[1], .x[2]))
conc_tbl <- map_dfr(conc, "table")
message("\n=== 2. Direction agreement between datasets (chance = 50%) ===")
print(conc_tbl, n = Inf, width = Inf)
write_csv(conc_tbl, file.path(in_dir, "diagnostics_concordance.csv"))

p_scatter <- conc[[1]]$merged |>
  ggplot(aes(stat_a, stat_b)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "#2166ac") +
  labs(x = "t-statistic, GSE67528 neurons", y = "t-statistic, GSE124308 neurons",
       title = sprintf("Neuron ASD effects across datasets (Spearman rho = %.3f)",
                       conc_tbl$spearman_t[1])) +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "neuron_concordance.png"), p_scatter,
       width = 6.5, height = 6, dpi = 300, bg = "white")

# --- 3. The source-adjusted hits across analyses -------------------------------
hits <- res$GSE124308_both_adjSource |> filter(padj < 0.05) |> pull(Gene_ID)
hit_tbl <- imap_dfr(res, ~ .x |> filter(Gene_ID %in% hits) |>
                      transmute(Gene_ID, Gene_Symbol, analysis = .y,
                                log2FC = round(log2FoldChange, 2),
                                pvalue = signif(pvalue, 2))) |>
  pivot_wider(id_cols = c(Gene_ID, Gene_Symbol), names_from = analysis,
              values_from = c(log2FC, pvalue))
message("\n=== 3. Source-adjusted hits in every analysis ===")
print(hit_tbl, n = Inf, width = Inf)
write_csv(hit_tbl, file.path(in_dir, "diagnostics_adjSource_hits.csv"))
