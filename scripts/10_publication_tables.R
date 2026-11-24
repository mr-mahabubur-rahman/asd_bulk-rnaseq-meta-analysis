# ==============================================================================
# 10_publication_tables.R — PUBLICATION-READY TABLES (discovery paper)
# ------------------------------------------------------------------------------
# Main tables  -> results/publication/Main_Tables.docx (landscape, Word-editable)
#   Table 1  Discovery cohorts: design, donors, samples, model
#   Table 2  Direction replication (H1): primary + all sensitivity analyses
#            (includes a FOXJ1-adjusted replication, computed here from the
#            existing 06 output — exploratory)
#   Table 3  Gene-set results by cohort (fgsea NES + cameraPR FDR)
#   Table 4  Motile/primary cilia sets before and after FOXJ1 adjustment (H2)
# Supplementary tables -> results/publication/Supplementary_Tables.xlsx
#   README sheet + S1–S10
# Supplementary data   -> results/publication/Supplementary_Data_1_DE_results.xlsx
#   full per-gene results for every model + the neuron meta-analysis
# Every table is also written as CSV (results/publication/csv/).
#
# Reads existing result tables only. Needs: r-flextable r-officer r-openxlsx
#   pixi add r-flextable r-officer r-openxlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(flextable)
  library(officer)
  library(openxlsx)
})
out_dir <- "results/publication"
csv_dir <- file.path(out_dir, "csv")
for (d in c(out_dir, csv_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- formatting helpers ---------------------------------------------------------
sup_map <- c("-" = "⁻", "0" = "⁰", "1" = "¹", "2" = "²", "3" = "³",
             "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸",
             "9" = "⁹")
sup <- function(e) paste(sup_map[strsplit(as.character(e), "")[[1]]], collapse = "")
sci1 <- function(p) {
  e <- floor(log10(p)); m <- round(p / 10^e, 1)
  if (m >= 10) { m <- m / 10; e <- e + 1 }
  paste0(formatC(m, format = "f", digits = 1), " × 10", sup(e))
}
fmt_p <- function(p) vapply(p, function(x) {
  if (is.na(x)) return("—")
  if (x >= 0.01) return(formatC(x, format = "f", digits = 2))
  if (x >= 0.001) return(formatC(x, format = "f", digits = 3))
  sci1(x)
}, character(1))
fmt_num <- function(x, d = 2) ifelse(is.na(x), "—",
                                     sub("^-", "−", formatC(x, format = "f", digits = d)))
fmt_pct <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)

pretty <- c(GSE67528_Neuron = "GSE67528 neurons", GSE67528_NPC = "GSE67528 NPCs",
            GSE124308_both = "GSE124308 neurons", GSE124308_no101_102 = "GSE124308 neurons (without 101/102)",
            GSE142670_Neuron = "GSE142670 neurons (mixed)",
            GSE142670_Neuron_sorted = "GSE142670 neurons (FACS-sorted)",
            GSE142670_NPC = "GSE142670 NPCs", pooled_neuron = "Pooled neurons",
            pooled_npc = "Pooled NPCs")

ft_style <- function(ft, widths = NULL) {
  ft <- ft |> theme_booktabs() |> font(fontname = "Arial", part = "all") |>
    fontsize(size = 8.5, part = "body") |> fontsize(size = 8.5, part = "header") |>
    fontsize(size = 7.5, part = "footer") |> bold(part = "header") |>
    valign(valign = "top", part = "all") |> padding(padding = 2, part = "all")
  if (is.null(widths)) set_table_properties(ft, layout = "autofit", width = 1)
  else width(ft, width = widths)
}
add_notes <- function(ft, notes) add_footer_lines(ft, notes) |>
  fontsize(size = 7.5, part = "footer") |> font(fontname = "Arial", part = "footer") |>
  border(part = "footer", border = fp_border(width = 0))
save_csv <- function(df, name) write_csv(df, file.path(csv_dir, paste0(name, ".csv")))

donors <- read_csv("results/tables/DE_dream/summary.csv", show_col_types = FALSE) |>
  select(analysis, donors) |> deframe()

# ================================================================================
# TABLE 1 — discovery cohorts
# ================================================================================
message("Table 1 ...")
cohort_rows <- tribble(
  ~gse,        ~stage,          ~label,                 ~model,
  "GSE67528",  "NPC",           "GSE67528_NPC",         "~ condition + (1 | donor)",
  "GSE67528",  "Neuron",        "GSE67528_Neuron",      "~ condition + (1 | donor)",
  "GSE124308", NA,              "GSE124308_both",       "~ condition + time point + (1 | donor)",
  "GSE142670", "NPC",           "GSE142670_NPC",        "~ condition + (1 | donor)",
  "GSE142670", "Neuron",        "GSE142670_Neuron",     "~ condition + (1 | donor)",
  "GSE142670", "Neuron_sorted", "GSE142670_Neuron_sorted", "~ condition + (1 | donor)")
refs <- c(GSE67528 = "Marchetto et al. 2017 [1]", GSE124308 = "DeRosa et al. 2018 [2]",
          GSE142670 = "Griesi-Oliveira et al. 2021 [3]")
design <- c(GSE67528 = "2D; iPSC, NPC, neurons", GSE124308 = "2D cortical neurons; day 35 and 135",
            GSE142670 = "2D; NPC, mixed neuronal cultures, FACS-sorted neurons")
counts_src <- c(GSE67528 = "Series count matrix (GEO)", GSE124308 = "Series count matrix (GEO)",
                GSE142670 = "Authors' per-sample counts (GEO)")
stage_lab <- c(NPC = "NPC", Neuron = "Neurons (mixed culture)", Neuron_sorted = "Neurons (FACS-sorted)")

t1 <- pmap_dfr(cohort_rows, function(gse, stage, label, model) {
  cn <- names(read_csv(paste0("data/raw_counts_primary/", gse, "_raw_counts.csv"),
                       n_max = 1, show_col_types = FALSE))
  m <- read_csv(paste0("data/metadata/", gse, "_metadata_full.csv"), show_col_types = FALSE) |>
    mutate(donor = as.character(donor)) |> filter(sample %in% cn)
  if (!is.na(stage)) m <- filter(m, stage == !!stage)
  per_donor <- count(m, donor)$n
  stage_txt <- if (is.na(stage)) {
    tp <- table(factor(m$timepoint, levels = c("d35", "d135")))
    paste0("Neurons (", paste(names(tp), tp, sep = ": n = ", collapse = "; "), ")")
  } else stage_lab[[stage]]
  tibble(Cohort = gse, Reference = refs[[gse]], Design = design[[gse]],
         `Stage analysed` = stage_txt,
         `Donors, ASD / control` = paste0(n_distinct(m$donor[m$condition == "ASD"]), " / ",
                                          n_distinct(m$donor[m$condition == "Control"])),
         `Samples, ASD / control` = paste0(sum(m$condition == "ASD"), " / ",
                                           sum(m$condition == "Control")),
         `Samples per donor` = if (min(per_donor) == max(per_donor)) as.character(min(per_donor))
                               else paste0(min(per_donor), "–", max(per_donor)),
         `Count data` = counts_src[[gse]], `Model (dream)` = model)
})
save_csv(t1, "Table1_cohorts")
# show cohort-level columns once per cohort (blank, not merged: merged cells
# draw stray borders in some Word/LibreOffice versions)
t1_show <- t1 |> group_by(Cohort) |>
  mutate(across(c(Reference, Design, `Count data`), \(x) if_else(row_number() == 1, x, ""))) |>
  ungroup() |> mutate(Cohort = if_else(duplicated(Cohort), "", Cohort))
ft1 <- flextable(t1_show) |>
  ft_style() |>
  add_notes(c(
    "All analyses used Ensembl release 115 primary-assembly genes (chromosomes 1–22, X, Y, MT); genes on alternative haplotypes, patches and retired IDs were removed. Each stage was analysed separately; donor was a random effect (variancePartition dream, voom weights, TMM normalisation, filterByExpr).",
    "GSE124308: all ASD donors were PBMC-derived; three of five control donors were fibroblast-derived (see Table 2 and Supplementary Table S8 for sensitivity analyses).",
    "FACS-sorted GSE142670 neurons were used only in a sensitivity analysis. NPC, neural progenitor cell; PBMC, peripheral blood mononuclear cell."))

# ================================================================================
# TABLE 2 — replication (H1)
# ================================================================================
message("Table 2 ...")
stouffer_z <- function(lfc, p, w) {
  z <- sign(lfc) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
  sum(w * z) / sqrt(sum(w^2))
}
load_de <- function(a, dir) read_csv(file.path(dir, paste0(a, ".csv")), show_col_types = FALSE) |>
  filter(is.finite(stat), is.finite(pvalue)) |> select(Gene_ID, log2FoldChange, pvalue) |>
  mutate(dataset = a)
replicate_loo <- function(long) {
  sets <- unique(long$dataset)
  map_dfr(sets, function(held) {
    top <- long |> filter(dataset != held) |> group_by(Gene_ID) |>
      filter(n() == length(sets) - 1) |>
      summarise(z = stouffer_z(log2FoldChange, pvalue, sqrt(donors[dataset])), .groups = "drop") |>
      arrange(desc(abs(z))) |> head(500)
    h <- long |> filter(dataset == held) |> inner_join(top, by = "Gene_ID")
    tibble(held_out = held, genes = nrow(h), k = sum(sign(h$log2FoldChange) == sign(h$z)))
  })
}
foxj1_dir <- "results/tables/DE_dream_adjFOXJ1"
neuron_sets <- c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron")
rep_foxj1 <- if (all(file.exists(file.path(foxj1_dir, paste0(neuron_sets, ".csv"))))) {
  replicate_loo(map_dfr(neuron_sets, load_de, dir = foxj1_dir)) |>
    mutate(analysis = "FOXJ1-adjusted (exploratory)")
} else NULL
if (!is.null(rep_foxj1)) {
  message("  FOXJ1-adjusted replication:")
  print(rep_foxj1 |> mutate(share = round(k / genes, 3)))
}

f2 <- read_csv("results/tables/figures/Fig2_replication.csv", show_col_types = FALSE)
c1 <- read_csv("results/tables/sensitivity/check1_replication.csv", show_col_types = FALSE) |>
  filter(version == "GSE124308 without 101/102") |>
  transmute(held_out, genes, k = round(share * genes),
            analysis = "Excluding GSE124308 control donors 101 and 102")
t2_raw <- bind_rows(
  f2 |> filter(panel == "A", version == "Unadjusted", held_out != "GSE142670_Neuron_sorted") |>
    transmute(held_out, genes, k, analysis = "Primary: donor-aware model"),
  f2 |> filter(panel == "A", version == "Neuron-score adjusted") |>
    transmute(held_out, genes, k, analysis = "Adjusted for neuronal marker score"),
  if (!is.null(rep_foxj1)) select(rep_foxj1, held_out, genes, k, analysis),
  c1,
  f2 |> filter(held_out == "GSE142670_Neuron_sorted") |>
    transmute(held_out, genes, k, analysis = "FACS-sorted GSE142670 neurons replacing mixed cultures"),
  f2 |> filter(panel == "B", version == "Unadjusted") |>
    transmute(held_out, genes, k, analysis = "NPCs: donor-aware model"),
  f2 |> filter(panel == "B", version == "Neuron-score adjusted") |>
    transmute(held_out, genes, k, analysis = "NPCs: adjusted for neuronal marker score")) |>
  mutate(analysis = factor(analysis, levels = unique(analysis)),
         held_out = factor(held_out, levels = c("GSE67528_Neuron", "GSE124308_both",
                                                "GSE124308_no101_102", "GSE142670_Neuron",
                                                "GSE142670_Neuron_sorted", "GSE67528_NPC",
                                                "GSE142670_NPC"))) |>
  arrange(analysis, held_out) |>
  mutate(analysis = as.character(analysis), held_out = as.character(held_out)) |>
  rowwise() |>
  mutate(share = k / genes,
         lo = binom.test(k, genes)$conf.int[1], hi = binom.test(k, genes)$conf.int[2],
         p1 = binom.test(k, genes, 0.5, alternative = "greater")$p.value,
         p2 = binom.test(k, genes, 0.5)$p.value) |> ungroup()
save_csv(t2_raw, "Table2_replication")
t2 <- t2_raw |>
  transmute(Analysis = analysis,
            `Held-out cohort` = coalesce(pretty[held_out], held_out),
            `Genes tested` = genes,
            `Same direction, n (%)` = paste0(k, " (", fmt_pct(share), ")"),
            `95% CI (%)` = paste0(fmt_pct(lo), "–", fmt_pct(hi)),
            `One-sided p` = fmt_p(p1), `Two-sided p` = fmt_p(p2))
ft2 <- flextable(t2) |> merge_v(j = "Analysis") |> ft_style() |> fix_border_issues() |>
  hline(i = which(diff(as.integer(factor(t2$Analysis, levels = unique(t2$Analysis)))) != 0),
        border = fp_border(color = "grey70", width = 0.5)) |>
  add_notes(c(
    "For each held-out cohort, the 500 genes with the largest |weighted Stouffer z| (weights √donors) in the remaining cohort(s) were selected, and agreement in direction of change in the held-out cohort was counted. Fewer than 500 genes were testable where the held-out cohort did not express all selected genes. 95% CI, exact binomial (Clopper–Pearson). One-sided p tests agreement > 50% (H1 is directional); two-sided p is shown because the sorted-neuron agreement is below 50%.",
    "FOXJ1-adjusted: models with FOXJ1 expression as a covariate (exploratory). Adjusting for a single gene that tracks the signal of interest can over-correct, so agreement below 50% after adjustment should not be read as opposite biology. Sorted-neuron row: genes selected as the top 500 by Stouffer p in the other two cohorts (FDR < 0.05 rule not met); sorted neurons came from three ASD and three control donors.",
    "NPCs: two cohorts only, so each cohort’s top 500 genes were tested in the other."))

# ================================================================================
# TABLE 3 — gene sets by cohort
# ================================================================================
message("Table 3 ...")
f3 <- read_csv("results/tables/figures/Fig3_enrichment_fgsea_cameraPR.csv", show_col_types = FALSE)
col_order <- c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron", "pooled_neuron",
               "GSE67528_NPC", "GSE142670_NPC", "pooled_npc")
col_lab <- c(GSE67528_Neuron = "GSE67528", GSE124308_both = "GSE124308",
             GSE142670_Neuron = "GSE142670", pooled_neuron = "Pooled",
             GSE67528_NPC = "GSE67528 ", GSE142670_NPC = "GSE142670 ", pooled_npc = "Pooled ")
set_order <- unique(f3$label[order(match(f3$theme, unique(f3$theme)))])
t3_long <- f3 |>
  mutate(cell = paste0(fmt_num(NES), " (", fmt_p(camera_fdr), ")"),
         sig = !is.na(camera_fdr) & camera_fdr < 0.05,
         ranking = factor(ranking, levels = col_order),
         label = factor(label, levels = set_order))
save_csv(t3_long |> select(theme, label, ranking, NES, fgsea_p = pval, camera_dir, camera_p,
                           camera_fdr), "Table3_gene_sets_long")
t3 <- t3_long |> select(Theme = theme, `Gene set` = label, ranking, cell) |>
  arrange(ranking) |> pivot_wider(names_from = ranking, values_from = cell) |>
  arrange(`Gene set`) |> mutate(across(everything(), as.character))
sigm <- t3_long |> select(label, ranking, sig) |> arrange(ranking) |>
  pivot_wider(names_from = ranking, values_from = sig) |> arrange(label)
names(t3) <- c("Theme", "Gene set", col_lab[names(t3)[-(1:2)]])
t3$Theme <- sub("†", "", t3$Theme)
ft3 <- flextable(t3) |>
  add_header_row(values = c("", "", "Neurons", "NPCs"), colwidths = c(1, 1, 4, 3)) |>
  merge_v(j = "Theme") |> ft_style() |> fix_border_issues() |>
  align(i = 1, align = "center", part = "header")
for (jj in seq_along(col_order)) {
  rows <- which(sigm[[col_order[jj]]] %in% TRUE)
  if (length(rows)) ft3 <- bold(ft3, i = rows, j = jj + 2, part = "body")
}
ft3 <- add_notes(ft3, c(
  "Cells show the fgsea normalised enrichment score (NES; positive, higher in ASD) and, in parentheses, the cameraPR false discovery rate (inter-gene correlation 0.05; Benjamini–Hochberg across the 12 sets within each column). Bold: cameraPR FDR < 0.05. Genes were ranked by the moderated t of the ASD coefficient (per cohort) or by the weighted Stouffer z (pooled columns; not independent of the cohort columns).",
  "Interferon sets were pre-specified because of the earlier pilot analysis. The GSE124308 interferon-α result was robust to adjustment for iPSC source and to excluding control donors 101 and 102 (Supplementary Table S8)."))

# ================================================================================
# TABLE 4 — cilia before/after FOXJ1 adjustment
# ================================================================================
message("Table 4 ...")
f4 <- read_csv("results/tables/figures/Fig4_cilia_adjustment.csv", show_col_types = FALSE)
t4_long <- f4 |>
  mutate(version = recode(version, `+ motile-cilia score (circular check)` = "+ motile-cilia score†"),
         cell = paste0(fmt_p(p_lower_one_sided), " (", fmt_num(NES), ")"))
save_csv(t4_long, "Table4_cilia_adjustment_long")
t4 <- t4_long |>
  mutate(Cohort = coalesce(pretty[analysis], analysis),
         set = sub(" \\(.*", "", set)) |>
  mutate(set = factor(set, levels = c("Motile-cilia set", "Primary-cilia set")),
         analysis = factor(analysis, levels = c("GSE67528_Neuron", "GSE124308_both",
                                                "GSE142670_Neuron", "GSE142670_NPC"))) |>
  arrange(set, analysis) |>
  mutate(set = as.character(set)) |>
  select(`Gene set` = set, Cohort, version, cell) |>
  pivot_wider(names_from = version, values_from = cell)
ft4 <- flextable(t4) |> merge_v(j = "Gene set") |> ft_style() |> fix_border_issues() |>
  add_notes(c(
    "Cells: one-sided cameraPR p for the set being lower in ASD (inter-gene correlation 0.05) and, in parentheses, the fgsea NES. Covariates were added to the donor-aware model; GSE67528 NPCs were not included in this analysis.",
    "†The motile-cilia score is built from the genes of the motile-cilia set, so adjusting for it removes the set’s shared variation by construction; it shows that the signal is a coordinated, programme-level shift and is not independent evidence.",
    "Curated sets: motile-cilia set (49 genes; axonemal dyneins, radial spoke, central pair and FOXJ1/MCIDAS programme); primary-cilia set (53 genes; intraflagellar transport, BBSome, transition zone, ciliary signalling)."))

# --- write the Word document ------------------------------------------------------
title_par <- function(txt) fpar(ftext(txt, fp_text(bold = TRUE, font.family = "Arial", font.size = 10)))
doc <- read_docx() |>
  body_set_default_section(prop_section(page_size = page_size(orient = "landscape"),
                                        page_margins = page_mar(top = 0.8, bottom = 0.8,
                                                                left = 0.8, right = 0.8)))
add_tab <- function(doc, title, ft, last = FALSE) {
  doc <- doc |> body_add_fpar(title_par(title)) |> body_add_flextable(ft, align = "left")
  if (!last) doc <- body_add_break(doc)
  doc
}
doc <- doc |>
  add_tab("Table 1. Discovery cohorts and analysis design", ft1) |>
  add_tab("Table 2. Direction replication of ASD-associated expression changes across cohorts (H1) and sensitivity analyses", ft2) |>
  add_tab("Table 3. Selected gene sets by cohort and developmental stage", ft3) |>
  add_tab("Table 4. Cilia gene sets before and after adjustment for FOXJ1 expression (H2)", ft4, last = TRUE)
print(doc, target = file.path(out_dir, "Main_Tables.docx"))
message("  wrote Main_Tables.docx")

# ================================================================================
# SUPPLEMENTARY TABLES (xlsx)
# ================================================================================
message("Supplementary tables ...")
rd <- function(path) if (file.exists(path)) read_csv(path, show_col_types = FALSE) else
  tibble(note = paste("file not found:", path))
# round only non-integer columns (signif() on counts such as 14574 would give 14570)
round_df <- function(df, d = 4) df |>
  mutate(across(where(is.double), \(x) ifelse(!is.na(x) & x == round(x), x, signif(x, d))))

# readable column headers for the supplement
hdr <- c(analysis = "Analysis", dataset = "Analysis", samples = "Samples", donors = "Donors",
         genes_tested = "Genes tested", padj_lt_0.05 = "Genes FDR < 0.05",
         padj_lt_0.05_lfc_gt_1 = "Genes FDR < 0.05 and |log2FC| > 1", min_padj = "Smallest FDR",
         min_fdr = "Smallest FDR", pi1 = "π1 (share non-null)", p_lt_0.01 = "Genes p < 0.01",
         expected_0.01 = "Expected p < 0.01 (null)", datasets = "Cohorts",
         fdr_lt_0.05 = "Genes FDR < 0.05 (random effects)", robust = "Robust genes",
         stouffer_fdr_lt_0.05 = "Genes FDR < 0.05 (Stouffer)", Gene_ID = "Ensembl gene ID",
         Gene_Symbol = "Gene", symbol = "Gene", gene = "Gene", k = "Cohorts",
         same_sign = "Same direction in all cohorts", est = "Pooled log2FC", se = "SE",
         ci_lb = "95% CI lower", ci_ub = "95% CI upper", p = "p", fdr = "FDR", tau2 = "τ²",
         I2 = "I² (%)", st_z = "Stouffer z", st_p = "Stouffer p", stouffer_fdr = "Stouffer FDR",
         score = "Marker score", n_asd = "ASD donors", n_ctrl = "Control donors",
         diff_asd_minus_ctrl = "ASD − control (SD)", t_p = "t-test p", wilcox_p = "Wilcoxon p",
         label = "Gene set", theme = "Theme", ranking = "Ranking", NES = "NES (fgsea)",
         pval = "fgsea p", camera_dir = "cameraPR direction", camera_p = "cameraPR p",
         camera_fdr = "cameraPR FDR", inter_gene_cor = "Inter-gene correlation", set = "Gene set",
         n_genes = "Genes in set", direction = "Direction", mean_logCPM = "Mean log2 CPM",
         asd_minus_ctrl = "ASD − control (log2 CPM)", p_donor_t = "Donor-level t-test p",
         intersection = "Intersection", observed = "Observed", expected = "Expected",
         ratio = "Observed / expected", perm_p_greater = "Permutation p", universe = "Genes in universe",
         pooled_z = "Pooled z", cilia = "Cilia-annotated", stage = "Culture",
         mean_ctrl = "Control mean (SD)", mean_asd = "ASD mean (SD)",
         mean_shift = "Sorted − mixed (SD)", p_paired = "Paired t-test p", pair = "Comparison",
         spearman_all = "Spearman ρ, all genes", spearman_top500 = "Spearman ρ, top-500 genes",
         asd_dir_other = "Direction in other cohorts", genes = "Genes",
         mean_sorting_t = "Mean sorting t", share_depleted_by_sorting = "Share depleted by sorting",
         share_same_dir_mixed = "Same direction, mixed cultures",
         share_same_dir_sorted = "Same direction, sorted neurons", metric = "Metric", value = "Value")
nice <- function(df) { names(df) <- coalesce(unname(hdr[names(df)]), names(df)); df }

c1_cam <- rd("results/tables/sensitivity/check1_cameraPR.csv")
supp <- list(
  S1_DE_summary = list(
    "Per-model differential expression summary: samples, donors, genes tested, FDR hits and π1 (share of non-null genes). Primary models: GSE67528_NPC, GSE67528_Neuron, GSE124308_both, GSE142670_NPC, GSE142670_Neuron; the others are sensitivity or per-time-point models.",
    rd("results/tables/DE_dream/summary.csv") |>
      left_join(rd("results/tables/DE_dream/diagnostics_pi1.csv") |>
                  select(any_of(c("analysis", "pi1", "p_lt_0.01", "expected_0.01"))),
                by = "analysis")),
  S2_meta_summary = list(
    "Random-effects meta-analysis of neurons (REML, Hartung–Knapp ad hoc): primary (mixed GSE142670 neurons) and sensitivity (sorted neurons). With three cohorts the test has two degrees of freedom and is very conservative.",
    rd("results/tables/meta/meta_summary.csv")),
  S3_meta_top100 = list(
    "Top 100 genes of the primary neuron meta-analysis by random-effects p (none reaches FDR < 0.05).",
    rd("results/tables/meta/neuron_primary_REM.csv") |> arrange(p) |> head(100) |>
      select(any_of(c("Gene_Symbol", "Gene_ID", "k", "same_sign", "est", "se", "ci_lb", "ci_ub",
                      "p", "fdr", "tau2", "I2", "st_z", "st_p", "stouffer_fdr")))),
  S4_composition = list(
    "Cell-composition marker scores, ASD vs control on donor means (difference in SD units; t-test and Wilcoxon).",
    rd("results/tables/composition/composition_tests_donor_level.csv")),
  S5_gene_sets_all = list(
    "fgsea and cameraPR results for the 12 gene sets of Table 3, all rankings (pooled rankings use the weighted Stouffer z).",
    rd("results/tables/figures/Fig3_enrichment_fgsea_cameraPR.csv") |>
      select(any_of(c("theme", "label", "ranking", "NES", "pval", "camera_dir", "camera_p",
                      "camera_fdr")))),
  S6_cilia_cameraPR = list(
    "cameraPR per cohort for GO cilia sets, curated motile and primary cilia sets and comparison sets, at inter-gene correlation 0.01 and 0.05.",
    rd("results/tables/cilia/cameraPR_per_dataset.csv")),
  S7_ependymal_markers = list(
    "Ependymal / choroid-plexus marker expression (FOXJ1, TTR, AQP1, CLIC6, KCNJ13, S100B) by cohort and diagnosis.",
    rd("results/tables/cilia/ependymal_markers.csv")),
  S8_sensitivity_GSE124308 = list(
    "GSE124308 without control donors 101 and 102 (fibroblast-derived; separate on PC2, Fig. S1B).",
    list(
      "A. Gene-level agreement between the full model and the model without 101/102" =
        rd("results/tables/sensitivity/check1_gene_level.csv") |>
        pivot_longer(everything(), names_to = "metric", values_to = "value") |>
        mutate(metric = recode(metric, genes_both = "Genes in both models",
                               spearman_t = "Spearman ρ of moderated t",
                               same_sign_share = "Share with the same direction",
                               p_lt_0.05_full = "Genes p < 0.05, full model",
                               p_lt_0.05_reduced = "Genes p < 0.05, without 101/102",
                               fdr_lt_0.05_reduced = "Genes FDR < 0.05, without 101/102")),
      "B. Leave-one-cohort-out replication with the reduced GSE124308 model" =
        rd("results/tables/sensitivity/check1_replication.csv"),
      "C. cameraPR (inter-gene correlation 0.05): direction in ASD and p" =
        c1_cam |> mutate(cell = paste0(direction, ", p = ", signif(p, 2))) |>
        select(set, version, cell) |> pivot_wider(names_from = version, values_from = cell))),
  S9_sorted_vs_mixed = list(
    "GSE142670 FACS-sorted vs mixed neuronal cultures (sorted neurons: 3 ASD and 3 control donors).",
    list(
      "A. Composition scores by diagnosis (SD units across all neuron samples; donor means)" =
        rd("results/tables/sensitivity/check2_composition_asd_vs_ctrl.csv"),
      "B. What sorting changes: sorted minus mixed, same donors (n = 6)" =
        rd("results/tables/sensitivity/check2_composition_sorted_minus_mixed.csv"),
      "C. Correlations between effects (sorting effect: sorted vs mixed, donor random effect)" =
        rd("results/tables/sensitivity/check2_correlations.csv"),
      "D. The 500 replication genes, by their direction in GSE67528 and GSE124308 neurons" =
        rd("results/tables/sensitivity/check2_top500_by_direction.csv"))),
  S10_top_overlaps = list(
    "Same-direction overlap of the top 500 genes between neuron cohorts, observed vs expected (2,000 permutations).",
    rd("results/tables/figures/FigS3_upset_observed_expected.csv")),
  S11_top40_genes = list(
    "The 40 genes with the strongest pooled neuron signal (same direction in all three neuron cohorts): log2FC and p per cohort. No gene reaches FDR < 0.05.",
    rd("results/tables/figures/FigS2_heatmap_top40.csv") |>
      select(any_of(c("symbol", "Gene_ID", "pooled_z", "cilia", "dataset", "log2FoldChange", "pvalue"))) |>
      pivot_wider(names_from = dataset, values_from = c(log2FoldChange, pvalue),
                  names_glue = "{.value} | {dataset}") |>
      rename_with(\(x) sub("^log2FoldChange \\| ", "log2FC, ", sub("^pvalue \\| ", "p, ", x))) |>
      arrange(desc(pooled_z)))
)

hs <- createStyle(textDecoration = "bold", fgFill = "#E7E6E6", border = "bottom", wrapText = TRUE)
bs <- createStyle(textDecoration = "bold")
wb <- createWorkbook()
addWorksheet(wb, "README")
readme <- tibble(Sheet = names(supp), Description = map_chr(supp, 1))
writeData(wb, "README", "Supplementary Tables — donor-aware re-analysis of iPSC-derived neural transcriptomes in idiopathic ASD (discovery cohorts GSE67528, GSE124308, GSE142670)", startRow = 1)
addStyle(wb, "README", bs, rows = 1, cols = 1)
writeData(wb, "README", readme, startRow = 3, headerStyle = hs)
setColWidths(wb, "README", cols = 1:2, widths = c(26, 120))
for (nm in names(supp)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, supp[[nm]][[1]], startRow = 1)
  addStyle(wb, nm, createStyle(textDecoration = "italic"), rows = 1, cols = 1)
  content <- supp[[nm]][[2]]
  if (is.data.frame(content)) {
    df <- round_df(content)
    save_csv(df, nm)
    df <- nice(df)
    writeData(wb, nm, df, startRow = 3, headerStyle = hs)
    freezePane(wb, nm, firstActiveRow = 4)
    setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
  } else {                                   # several blocks stacked on one sheet
    row <- 3; maxc <- 1
    for (bt in names(content)) {
      df <- round_df(content[[bt]])
      save_csv(df, paste0(nm, "_", substr(bt, 1, 1)))
      df <- nice(df)
      writeData(wb, nm, bt, startRow = row); addStyle(wb, nm, bs, rows = row, cols = 1)
      writeData(wb, nm, df, startRow = row + 1, headerStyle = hs)
      row <- row + nrow(df) + 3; maxc <- max(maxc, ncol(df))
    }
    setColWidths(wb, nm, cols = seq_len(maxc), widths = "auto")
  }
}
saveWorkbook(wb, file.path(out_dir, "Supplementary_Tables.xlsx"), overwrite = TRUE)
message("  wrote Supplementary_Tables.xlsx")

# ================================================================================
# SUPPLEMENTARY DATA 1 — full per-gene results
# ================================================================================
message("Supplementary Data 1 ...")
de_files <- list.files("results/tables/DE_dream", pattern = "^GSE.*\\.csv$", full.names = TRUE)
wb2 <- createWorkbook()
addWorksheet(wb2, "README")
writeData(wb2, "README", tibble(
  Sheet = c(tools::file_path_sans_ext(basename(de_files)), "Meta_neurons_REM"),
  Description = c(rep("Per-gene ASD vs control results (dream; log2FoldChange, lfcSE, moderated t, p, BH FDR).",
                      length(de_files)),
                  "Per-gene random-effects meta-analysis of the three neuron cohorts.")),
  headerStyle = hs)
writeData(wb2, "README", tibble(
  Column = c("Gene_ID", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "AveExpr",
             "Gene_Symbol / Gene_Biotype", "Meta: est, se, ci_lb, ci_ub, p, fdr",
             "Meta: tau2, I2, Q", "Meta: st_z, st_p, stouffer_fdr", "Meta: max_p_leave_one_out, robust"),
  Definition = c("Ensembl gene ID (release 115, primary assembly)",
                 "ASD vs control log2 fold change (dream, donor random effect)",
                 "Standard error implied by the moderated t (|log2FC / t|)",
                 "Moderated t statistic", "p value", "Benjamini–Hochberg FDR",
                 "Average log2 CPM", "Ensembl gene name and biotype",
                 "Random-effects pooled log2FC, SE, 95% CI, p and FDR (REML, Hartung–Knapp ad hoc)",
                 "Between-cohort variance, I² (%), Cochran's Q",
                 "Weighted Stouffer z (weights √donors), p and FDR",
                 "Largest p when any one cohort is dropped (genes with FDR < 0.10 only); robust flag")),
  startRow = length(de_files) + 5, headerStyle = hs)
setColWidths(wb2, "README", cols = 1:2, widths = c(34, 100))
for (f in de_files) {
  nm <- substr(tools::file_path_sans_ext(basename(f)), 1, 31)
  addWorksheet(wb2, nm)
  df <- read_csv(f, show_col_types = FALSE) |> round_df(5)
  writeData(wb2, nm, df, headerStyle = hs); freezePane(wb2, nm, firstActiveRow = 2)
}
addWorksheet(wb2, "Meta_neurons_REM")
writeData(wb2, "Meta_neurons_REM", rd("results/tables/meta/neuron_primary_REM.csv") |> round_df(5),
          headerStyle = hs)
freezePane(wb2, "Meta_neurons_REM", firstActiveRow = 2)
saveWorkbook(wb2, file.path(out_dir, "Supplementary_Data_1_DE_results.xlsx"), overwrite = TRUE)
message("  wrote Supplementary_Data_1_DE_results.xlsx")
message("\nDone. Files in ", out_dir)
