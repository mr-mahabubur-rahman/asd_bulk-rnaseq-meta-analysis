# ==============================================================================
# 09_sensitivity_checks.R — TWO CHECKS BEFORE WRITING THE DISCOVERY PAPER
# ------------------------------------------------------------------------------
# CHECK 1  GSE124308 without control donors 101 and 102
#          (fibroblast-derived; separate from every other sample on PC2, Fig S1B;
#          all ASD donors are PBMC-derived). Refit the same dream model and ask:
#            a. how similar are the gene-level results to the full model?
#            b. does leave-one-cohort-out replication change?
#            c. do the pathway results change (interferon, MYC, translation, cilia)?
#
# CHECK 2  Why do FACS-sorted GSE142670 neurons go the OPPOSITE way (33%)?
#            a. which donors / how many samples have sorted vs mixed neurons
#            b. composition scores, sorted vs mixed, ASD vs control (donor level)
#            c. "sorting effect" = sorted vs mixed in the same donors (dream,
#               donor random effect). If the ASD signature in the other cohorts
#               is made of genes that sorting removes (non-neuronal cells),
#               sorting changes what an ASD-vs-control contrast can see.
#
# All exploratory: reported as sensitivity analyses in the discovery paper.
# Output: results/tables/sensitivity/, results/tables/DE_dream_sens/, logs via tee
# Runtime: two dream fits (several minutes each)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
  library(msigdbr)
})
set.seed(2026)
param   <- MulticoreParam(workers = 4)
tab_dir <- "results/tables/sensitivity"
de_dir  <- "results/tables/DE_dream_sens"
for (d in c(tab_dir, de_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

genes <- read_tsv("data/annotation/ensembl115_genes.tsv",
                  col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE)
sym_of <- setNames(genes$external_gene_name, genes$ensembl_gene_id)
donors <- read_csv("results/tables/DE_dream/summary.csv", show_col_types = FALSE) |>
  select(analysis, donors) |> deframe()

read_counts <- function(gse) {
  x <- read_csv(paste0("data/raw_counts_primary/", gse, "_raw_counts.csv"),
                show_col_types = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m
}
read_meta <- function(gse) {
  read_csv(paste0("data/metadata/", gse, "_metadata_full.csv"), show_col_types = FALSE) |>
    mutate(condition = factor(condition, levels = c("Control", "ASD")),
           donor = factor(as.character(donor)))
}
load_de <- function(a, dir = "results/tables/DE_dream") {
  read_csv(file.path(dir, paste0(a, ".csv")), show_col_types = FALSE) |>
    filter(is.finite(stat), is.finite(pvalue)) |> mutate(dataset = a)
}
fit_dream <- function(counts, meta, form, coef) {
  meta <- droplevels(as.data.frame(meta)); rownames(meta) <- meta$sample
  counts <- counts[, meta$sample, drop = FALSE]
  dge  <- DGEList(counts)
  grp  <- if ("condition" %in% all.vars(form)) meta$condition else meta$stage
  dge  <- calcNormFactors(dge[filterByExpr(dge, group = grp), , keep.lib.sizes = FALSE])
  vobj <- voomWithDreamWeights(dge, form, meta, BPPARAM = param)
  fit  <- eBayes(dream(vobj, form, meta, BPPARAM = param))
  tt   <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tibble(Gene_ID = rownames(tt), log2FoldChange = tt$logFC, lfcSE = abs(tt$logFC / tt$t),
         stat = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val,
         Gene_Symbol = unname(sym_of[rownames(tt)]))
}
stouffer_z <- function(lfc, p, w) {
  z <- sign(lfc) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
  sum(w * z) / sqrt(sum(w^2))
}
# same rule as 03_cell_composition.R: top 500 by |weighted Stouffer z| in the others
replicate_loo <- function(tables, version) {
  long <- bind_rows(tables)
  sets <- unique(long$dataset)
  map_dfr(sets, function(held) {
    top <- long |> filter(dataset != held) |>
      group_by(Gene_ID) |> filter(n() == length(sets) - 1) |>
      summarise(z = stouffer_z(log2FoldChange, pvalue, sqrt(donors[dataset])),
                .groups = "drop") |>
      arrange(desc(abs(z))) |> head(500)
    h <- long |> filter(dataset == held) |> inner_join(top, by = "Gene_ID")
    k <- sum(sign(h$log2FoldChange) == sign(h$z))
    tibble(version = version, held_out = held, genes = nrow(h), share = k / nrow(h),
           p_one_sided = binom.test(k, nrow(h), 0.5, alternative = "greater")$p.value,
           p_two_sided = binom.test(k, nrow(h), 0.5)$p.value)
  })
}

# gene sets for the pathway comparison
get_sets <- function(coll, sub = NULL) {
  new_api <- "collection" %in% names(formals(msigdbr))
  df <- if (new_api) {
    if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
    else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  } else {
    if (is.null(sub)) msigdbr(species = "Homo sapiens", category = coll)
    else msigdbr(species = "Homo sapiens", category = coll, subcategory = sub)
  }
  split(df$gene_symbol, df$gs_name) |> map(unique)
}
all_sets <- c(get_sets("H"), get_sets("C5", "GO:BP"), get_sets("C2", "CP:REACTOME"))
motile <- c("DNAH1","DNAH2","DNAH3","DNAH5","DNAH6","DNAH7","DNAH9","DNAH10","DNAH11",
            "DNAH12","DNAI1","DNAI2","DNAI3","DNAI4","DNAL1","RSPH1","RSPH3","RSPH4A",
            "RSPH9","SPAG6","SPAG16","SPAG17","HYDIN","CFAP43","CFAP44","CCDC39","CCDC40",
            "ODAD1","ODAD2","DNAAF1","DNAAF2","DNAAF4","ZMYND10","FOXJ1","MCIDAS","GMNC",
            "CDC20B","TEKT1","TEKT2","TPPP3","CAPS","SNTN","DYNLRB2","PACRG","RIBC1",
            "RIBC2","CFAP45","CFAP52","ENKUR")
test_sets <- c(all_sets[intersect(c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                                    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                    "HALLMARK_MYC_TARGETS_V1",
                                    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
                                    "GOBP_AXONEME_ASSEMBLY"), names(all_sets))],
               list(CURATED_MOTILE_CILIA = motile))
camera_one <- function(res, version) {
  st <- res |> filter(!is.na(Gene_Symbol), is.finite(stat)) |>
    group_by(Gene_Symbol) |> slice_max(abs(stat), n = 1, with_ties = FALSE) |>
    ungroup() |> select(Gene_Symbol, stat) |> deframe() |> sort(decreasing = TRUE)
  cameraPR(st, ids2indices(test_sets, names(st)), inter.gene.cor = 0.05, sort = FALSE) |>
    as_tibble(rownames = "set") |>
    transmute(version, set, direction = Direction, p = PValue)
}

# ================================================================================
# CHECK 1 — GSE124308 without donors 101 and 102
# ================================================================================
message("\n==================== CHECK 1: GSE124308 without 101, 102 ====================")
c2 <- read_counts("GSE124308")
m2 <- read_meta("GSE124308") |> filter(sample %in% colnames(c2)) |>
  mutate(timepoint = factor(timepoint))
drop_donors <- c("101", "102")
m2r <- m2 |> filter(!donor %in% drop_donors)
message("Donors kept by group and iPSC source:")
print(m2r |> distinct(donor, condition, ipsc_source) |> count(condition, ipsc_source))

red <- fit_dream(c2, m2r, ~ condition + timepoint + (1 | donor), "conditionASD")
write_csv(red, file.path(de_dir, "GSE124308_both_no101_102.csv"))
donors["GSE124308_no101_102"] <- n_distinct(m2r$donor)

full <- load_de("GSE124308_both")
cmp1 <- inner_join(full |> select(Gene_ID, t_full = stat, lfc_full = log2FoldChange),
                   red |> select(Gene_ID, t_red = stat, lfc_red = log2FoldChange),
                   by = "Gene_ID")
gene_level <- tibble(
  genes_both = nrow(cmp1),
  spearman_t = cor(cmp1$t_full, cmp1$t_red, method = "spearman"),
  same_sign_share = mean(sign(cmp1$lfc_full) == sign(cmp1$lfc_red)),
  p_lt_0.05_full = sum(full$pvalue < 0.05), p_lt_0.05_reduced = sum(red$pvalue < 0.05),
  fdr_lt_0.05_reduced = sum(red$padj < 0.05, na.rm = TRUE))
message("\n1a. Gene level, full vs without 101/102:")
print(gene_level, width = Inf)

neuron_full <- list(load_de("GSE67528_Neuron"), load_de("GSE124308_both"),
                    load_de("GSE142670_Neuron"))
neuron_red  <- list(load_de("GSE67528_Neuron"),
                    red |> mutate(dataset = "GSE124308_no101_102"),
                    load_de("GSE142670_Neuron"))
rep1 <- bind_rows(replicate_loo(neuron_full, "full GSE124308"),
                  replicate_loo(neuron_red, "GSE124308 without 101/102"))
message("\n1b. Leave-one-cohort-out replication (neurons, top 500):")
print(rep1 |> mutate(share = round(share, 3), across(starts_with("p_"), \(x) signif(x, 2))),
      width = Inf)

src_file <- "results/tables/DE_dream/GSE124308_both_adjSource.csv"
cam1 <- bind_rows(
  camera_one(full, "full"),
  if (file.exists(src_file)) camera_one(read_csv(src_file, show_col_types = FALSE),
                                        "adjusted for iPSC source"),
  camera_one(red, "without 101/102"))
message("\n1c. cameraPR (inter-gene cor 0.05), GSE124308:")
print(cam1 |> mutate(cell = paste0(direction, " p=", signif(p, 2))) |>
        select(set, version, cell) |> pivot_wider(names_from = version, values_from = cell),
      width = Inf)

write_csv(gene_level, file.path(tab_dir, "check1_gene_level.csv"))
write_csv(rep1, file.path(tab_dir, "check1_replication.csv"))
write_csv(cam1, file.path(tab_dir, "check1_cameraPR.csv"))

# ================================================================================
# CHECK 2 — sorted vs mixed GSE142670 neurons
# ================================================================================
message("\n==================== CHECK 2: GSE142670 sorted vs mixed neurons ====================")
c3 <- read_counts("GSE142670")
m3 <- read_meta("GSE142670") |>
  filter(sample %in% colnames(c3), stage %in% c("Neuron", "Neuron_sorted")) |>
  mutate(stage = factor(stage, levels = c("Neuron", "Neuron_sorted")))

message("2a. Samples per donor and stage:")
print(m3 |> count(condition, donor, stage) |>
        pivot_wider(names_from = stage, values_from = n, values_fill = 0), n = Inf)
paired <- m3 |> distinct(donor, stage) |> count(donor) |> filter(n == 2) |> pull(donor)
message("Donors with both mixed and sorted neurons: ", length(paired), " of ",
        n_distinct(m3$donor))

# 2b. composition scores on a common scale (z across all mixed + sorted samples)
markers <- list(
  neuron        = c("MAP2", "STMN2", "SYP", "SNAP25", "RBFOX3", "SYT1", "GAP43",
                    "DCX", "TUBB3", "ELAVL4", "INA", "NEFL", "NEFM"),
  progenitor    = c("NES", "SOX2", "PAX6", "VIM", "HES5", "FABP7", "HES1", "NOTCH1"),
  proliferation = c("MKI67", "TOP2A", "CENPF", "CCNB1", "BIRC5", "PCNA", "MCM2"),
  glia          = c("GFAP", "AQP4", "S100B", "SLC1A3", "ALDH1L1", "OLIG2"),
  non_neural    = c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "TAGLN", "DCN", "LUM", "POSTN"),
  motile_cilia  = motile,
  glutamatergic = c("SLC17A7", "SLC17A6", "GRIN2B", "TBR1", "NEUROD6"),
  gabaergic     = c("GAD1", "GAD2", "SLC32A1", "DLX5", "DLX6"))
lc  <- cpm(calcNormFactors(DGEList(c3[, m3$sample])), log = TRUE, prior.count = 2)
sym <- sym_of[rownames(lc)]
zr  <- t(scale(t(lc)))
sc  <- imap_dfc(markers, function(g, nm) {
  idx <- which(sym %in% g & rowMeans(lc) > 1)
  tibble(!!nm := if (length(idx) < 2) NA_real_ else colMeans(zr[idx, , drop = FALSE],
                                                               na.rm = TRUE))
})
scores <- bind_cols(m3 |> select(sample, donor, condition, stage), sc)
write_csv(scores, file.path(tab_dir, "check2_composition_scores.csv"))

dm <- scores |> group_by(stage, donor, condition) |>
  summarise(across(all_of(names(markers)), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(all_of(names(markers)), names_to = "score")
comp2 <- dm |> group_by(stage, score) |>
  summarise(mean_ctrl = mean(value[condition == "Control"]),
            mean_asd  = mean(value[condition == "ASD"]),
            diff_asd_minus_ctrl = mean_asd - mean_ctrl,
            t_p = tryCatch(t.test(value ~ condition)$p.value, error = function(e) NA_real_),
            .groups = "drop")
message("\n2b. Composition scores (SD units across all neuron samples; donor means):")
print(comp2 |> mutate(across(where(is.double), \(x) signif(x, 2))) |>
        arrange(score, stage), n = Inf, width = Inf)
sort_shift <- dm |> filter(donor %in% paired) |>
  pivot_wider(names_from = stage, values_from = value) |>
  mutate(sorted_minus_mixed = Neuron_sorted - Neuron) |>
  group_by(score) |>
  summarise(mean_shift = mean(sorted_minus_mixed),
            p_paired = tryCatch(t.test(sorted_minus_mixed)$p.value, error = function(e) NA_real_),
            .groups = "drop")
message("Sorted minus mixed, same donors (what sorting changes):")
print(sort_shift |> mutate(across(where(is.double), \(x) signif(x, 2))), width = Inf)
write_csv(comp2, file.path(tab_dir, "check2_composition_asd_vs_ctrl.csv"))
write_csv(sort_shift, file.path(tab_dir, "check2_composition_sorted_minus_mixed.csv"))

# 2c. sorting effect (sorted vs mixed, donor random effect, adjusted for diagnosis)
message("\n2c. Fitting sorted vs mixed (~ stage + condition + (1|donor)) ...")
sort_eff <- fit_dream(c3, m3, ~ stage + condition + (1 | donor), "stageNeuron_sorted")
write_csv(sort_eff, file.path(de_dir, "GSE142670_sorted_vs_mixed.csv"))

others <- bind_rows(load_de("GSE67528_Neuron"), load_de("GSE124308_both")) |>
  group_by(Gene_ID) |> filter(n() == 2) |>
  summarise(z_other = stouffer_z(log2FoldChange, pvalue, sqrt(donors[dataset])),
            .groups = "drop")
mixed  <- load_de("GSE142670_Neuron")  |> select(Gene_ID, t_mixed = stat)
sorted <- load_de("GSE142670_Neuron_sorted") |> select(Gene_ID, t_sorted = stat)
j <- others |> inner_join(sort_eff |> select(Gene_ID, t_sorting = stat), by = "Gene_ID") |>
  inner_join(mixed, by = "Gene_ID") |> inner_join(sorted, by = "Gene_ID")
top <- j |> arrange(desc(abs(z_other))) |> head(500) |>
  mutate(asd_dir_other = if_else(z_other > 0, "higher in ASD (other cohorts)",
                                 "lower in ASD (other cohorts)"))

cors <- tibble(
  pair = c("ASD effect: mixed vs sorted (GSE142670)",
           "ASD effect: other cohorts vs mixed",
           "ASD effect: other cohorts vs sorted",
           "sorting effect vs ASD effect in other cohorts",
           "sorting effect vs ASD effect in sorted"),
  spearman_all = c(cor(j$t_mixed, j$t_sorted, method = "spearman"),
                   cor(j$z_other, j$t_mixed, method = "spearman"),
                   cor(j$z_other, j$t_sorted, method = "spearman"),
                   cor(j$t_sorting, j$z_other, method = "spearman"),
                   cor(j$t_sorting, j$t_sorted, method = "spearman")),
  spearman_top500 = c(cor(top$t_mixed, top$t_sorted, method = "spearman"),
                      cor(top$z_other, top$t_mixed, method = "spearman"),
                      cor(top$z_other, top$t_sorted, method = "spearman"),
                      cor(top$t_sorting, top$z_other, method = "spearman"),
                      cor(top$t_sorting, top$t_sorted, method = "spearman")))
message("\n2c. Correlations (all genes; and the 500 genes used for replication):")
print(cors |> mutate(across(where(is.double), \(x) round(x, 3))), width = Inf)

by_dir <- top |> group_by(asd_dir_other) |>
  summarise(genes = n(),
            mean_sorting_t = mean(t_sorting),
            share_depleted_by_sorting = mean(t_sorting < 0),
            share_same_dir_mixed  = mean(sign(t_mixed)  == sign(z_other)),
            share_same_dir_sorted = mean(sign(t_sorted) == sign(z_other)),
            .groups = "drop")
message("Top-500 replication genes, split by their direction in the other cohorts:")
print(by_dir |> mutate(across(where(is.double), \(x) round(x, 3))), width = Inf)
message("(sorting t < 0 = gene is lower after sorting, i.e. enriched in non-neuronal cells)")

write_csv(cors, file.path(tab_dir, "check2_correlations.csv"))
write_csv(by_dir, file.path(tab_dir, "check2_top500_by_direction.csv"))
message("\nDone. Tables in ", tab_dir)
