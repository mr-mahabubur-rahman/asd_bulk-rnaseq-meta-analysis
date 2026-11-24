# ==============================================================================
# 08_main_figures.R — THE TWO MAIN FIGURES THAT CARRY H1 AND H2
# ------------------------------------------------------------------------------
#  Fig 2  Replication (H1): leave-one-cohort-out direction agreement of the
#         top-500 genes, unadjusted vs neuron-adjusted, with exact 95% CI and
#         ONE-SIDED binomial p (H1 is directional: agreement > 50%).
#         A: neurons (+ sorted-neuron sensitivity bar)   B: NPC (2 cohorts)
#  Fig 4  Cilia (H2): cameraPR for the curated motile and primary cilia sets,
#         before and after adjusting for FOXJ1 or the motile-cilia score.
#
# Reads existing tables only (03, 02, 06); nothing is refitted. Takes seconds.
# Output: results/figures/submission/, results/tables/figures/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
fig_dir <- "results/figures/submission"
tab_dir <- "results/tables/figures"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = pdf_dev)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300,
         bg = "white")
  message("  saved ", name, " (.pdf, .png)")
}
pretty <- c(GSE67528_Neuron = "GSE67528\nneurons", GSE124308_both = "GSE124308\nneurons",
            GSE142670_Neuron = "GSE142670\nneurons",
            GSE142670_Neuron_sorted = "GSE142670\nsorted neurons‡",
            GSE67528_NPC = "GSE67528\nNPC", GSE142670_NPC = "GSE142670\nNPC")
fmt_p <- function(p) if_else(p < 0.001, formatC(p, format = "e", digits = 1),
                             formatC(p, format = "f", digits = 3))

# ================================================================================
# Fig 2 — replication (H1)
# ================================================================================
message("\n== Fig 2: replication ==")
add_stats <- function(df) {
  df |> mutate(k = round(share * genes)) |>
    rowwise() |>
    mutate(lo = binom.test(k, genes)$conf.int[1],
           hi = binom.test(k, genes)$conf.int[2],
           p_one_sided = binom.test(k, genes, 0.5, alternative = "greater")$p.value,
           p_two_sided = binom.test(k, genes, 0.5)$p.value) |>
    ungroup()
}
rep_n <- read_csv("results/tables/composition/replication_unadj_vs_adjNeuron.csv",
                  show_col_types = FALSE)
sorted_file <- "results/tables/meta/neuron_sens_sorted_replication.csv"
if (file.exists(sorted_file)) {
  rep_n <- bind_rows(rep_n,
    read_csv(sorted_file, show_col_types = FALSE) |>
      filter(held_out == "GSE142670_Neuron_sorted") |>
      transmute(version = "unadjusted", held_out, genes, share))
}
rep_npc <- read_csv("results/tables/composition/npc_agreement_unadj_vs_adjNeuron.csv",
                    show_col_types = FALSE)
rep <- bind_rows(mutate(rep_n, panel = "A"), mutate(rep_npc, panel = "B")) |>
  add_stats() |>
  mutate(version = factor(recode(version, unadjusted = "Unadjusted",
                                 `neuron-adjusted` = "Neuron-score adjusted"),
                          levels = c("Unadjusted", "Neuron-score adjusted")),
         cohort = factor(pretty[held_out], levels = pretty),
         lab = paste0(round(100 * share), "%")) |>
  # empty partner bar for the sorted-neuron sensitivity, so every bar has the same width
  complete(nesting(panel, held_out, cohort), version)
write_csv(rep |> filter(!is.na(share)) |>
            select(panel, held_out, version, genes, k, share, lo, hi, p_one_sided, p_two_sided),
          file.path(tab_dir, "Fig2_replication.csv"))
print(rep |> filter(!is.na(share)) |>
        select(panel, held_out, version, genes, share, p_one_sided, p_two_sided) |>
        mutate(share = round(share, 3), across(starts_with("p_"), \(x) signif(x, 2))), n = Inf)

rep_plot <- function(d, xlab) {
  ggplot(d, aes(cohort, 100 * share, fill = version)) +
    geom_hline(yintercept = 50, linetype = 2, colour = "grey40") +
    geom_col(position = position_dodge(0.8), width = 0.75, na.rm = TRUE) +
    geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), position = position_dodge(0.8),
                  width = 0.25, linewidth = 0.4, na.rm = TRUE) +
    geom_text(aes(y = 100 * hi + 1.5, label = lab), position = position_dodge(0.8),
              vjust = 0, size = 2.5, na.rm = TRUE) +
    scale_fill_manual(values = c(Unadjusted = "grey55", `Neuron-score adjusted` = "#4d9221"),
                      name = NULL, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = xlab, y = "Same direction as the other cohort(s) (%)") +
    theme_classic(base_size = 9) + theme(legend.position = "top")
}
p2 <- rep_plot(filter(rep, panel == "A"), "Held-out neuron cohort") +
  rep_plot(filter(rep, panel == "B"), "Held-out NPC cohort") +
  plot_layout(widths = c(4, 2), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "top")
save_fig(p2, "Fig2_replication", 8, 4.2)

# ================================================================================
# Fig 4 — cilia before and after FOXJ1 / motile-score adjustment (H2)
# ================================================================================
message("\n== Fig 4: cilia adjustment ==")
cmp <- read_csv("results/tables/cilia/cilia_adjusted_comparison.csv", show_col_types = FALSE) |>
  filter(set %in% c("CURATED_MOTILE_CILIA", "CURATED_PRIMARY_CILIA")) |>
  mutate(signed = if_else(camera_dir == "Down", -1, 1) * -log10(camera_p),
         p_lower_one_sided = if_else(camera_dir == "Down", camera_p / 2, 1 - camera_p / 2),
         set = factor(recode(set, CURATED_MOTILE_CILIA = "Motile-cilia set",
                             CURATED_PRIMARY_CILIA = "Primary-cilia set (comparison)"),
                      levels = c("Motile-cilia set",
                                 "Primary-cilia set (comparison)")),
         version = factor(recode(version, unadjusted = "Unadjusted", adj_foxj1 = "+ FOXJ1",
                                 adj_motile = "+ motile-cilia score (circular check)"),
                          levels = c("Unadjusted", "+ FOXJ1",
                                     "+ motile-cilia score (circular check)")),
         cohort = factor(pretty[analysis], levels = pretty))
write_csv(cmp |> select(set, analysis, version, camera_dir, camera_p, p_lower_one_sided, NES),
          file.path(tab_dir, "Fig4_cilia_adjustment.csv"))
message("Motile-cilia set, one-sided p (lower in ASD):")
print(cmp |> filter(str_starts(set, "Motile")) |>
        select(analysis, version, camera_dir, p_lower_one_sided, NES) |>
        mutate(p_lower_one_sided = signif(p_lower_one_sided, 2), NES = round(NES, 2)) |>
        pivot_wider(names_from = version, values_from = c(p_lower_one_sided, NES)),
      width = Inf)

lim <- max(abs(cmp$signed), 2) * 1.1
p4 <- ggplot(cmp, aes(cohort, signed, fill = version)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(0.05), ymax = -log10(0.05),
           fill = "grey93") +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_wrap(~ set, ncol = 1) +
  scale_fill_manual(values = c(Unadjusted = "grey55", `+ FOXJ1` = "#2c7fb8",
                               `+ motile-cilia score (circular check)` = "#7fcdbb"),
                    name = NULL) +
  scale_y_continuous(limits = c(-lim, lim)) +
  labs(x = NULL,
       y = expression("Signed " * -log[10] * italic(p) * " (cameraPR; < 0 = lower in ASD)")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0))
save_fig(p4, "Fig4_cilia_adjustment", 6, 5.5)

message("\nDone. Figures in ", fig_dir)
