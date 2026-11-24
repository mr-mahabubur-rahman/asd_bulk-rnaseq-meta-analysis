# ==============================================================================
# 07b_figures_revised.R — SUBMISSION VERSIONS OF THE DESCRIPTIVE FIGURES
# ------------------------------------------------------------------------------
# Replaces the PNGs from 07_figures.R (07 stays in the repo as the first draft).
# No titles/subtitles on figures: that text belongs in the figure legends.
# Every figure is saved as PDF (vector) and PNG (300 dpi).
#
#  Fig 3      Enrichment dot plot, neurons | NPC, cameraPR FDR < 0.05 outlined
#  Fig S1     PCA, five panels in one figure (+ PC1/PC2 vs proliferation score)
#  Fig S2     Heatmap of the top-40 pooled genes, NA shown grey, cilia genes marked
#  Fig S3     UpSet with observed vs expected for every intersection
#  Fig S4     STRING network (exploratory), IFT/transition-zone relabelled
#  Fig S5     Volcano of the pooled neuron random-effects meta-analysis
#
# Needs: r-patchwork  (pixi add r-patchwork)
# Output: results/figures/submission/, results/tables/figures/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ggrepel)
  library(UpSetR)
  library(igraph)
  library(ggraph)
  library(httr)
  library(patchwork)
})
set.seed(2026)
fig_dir <- "results/figures/submission"
tab_dir <- "results/tables/figures"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# cairo_pdf keeps Greek letters and symbols (α, ∩, ≤); plain pdf() drops them
pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = pdf_dev)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300,
         bg = "white")
  message("  saved ", name, " (.pdf, .png)")
}
col_cond <- c(Control = "#2166ac", ASD = "#b2182b")

genes <- read_tsv("data/annotation/ensembl115_genes.tsv",
                  col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE)
sym_of <- setNames(genes$external_gene_name, genes$ensembl_gene_id)
donors <- read_csv("results/tables/DE_dream/summary.csv", show_col_types = FALSE) |>
  select(analysis, donors) |> deframe()

neuron <- c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron")
npc    <- c("GSE67528_NPC", "GSE142670_NPC")
pretty <- c(GSE67528_Neuron = "GSE67528 neurons", GSE124308_both = "GSE124308 neurons",
            GSE142670_Neuron = "GSE142670 neurons", GSE67528_NPC = "GSE67528 NPC",
            GSE142670_NPC = "GSE142670 NPC", pooled_neuron = "Pooled neurons",
            pooled_npc = "Pooled NPC")
load_de <- function(a) read_csv(file.path("results/tables/DE_dream", paste0(a, ".csv")),
                                show_col_types = FALSE) |>
  filter(is.finite(stat), is.finite(pvalue)) |> mutate(dataset = a)

# --- gene sets (shared by several figures) ----------------------------------------
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
cilia_any <- unique(c(motile, all_sets$GOBP_CILIUM_ORGANIZATION,
                      all_sets$GOBP_CILIUM_MOVEMENT, all_sets$GOBP_AXONEME_ASSEMBLY))

# ================================================================================
# Fig 3 — enrichment dot plot
# ================================================================================
message("\n== Fig 3: enrichment dot plot ==")
themes <- tribble(
  ~set,                                         ~label,                                ~theme,
  "GOBP_CILIUM_ORGANIZATION",                   "Cilium organization (GO:BP)",         "Cilia",
  "GOBP_AXONEME_ASSEMBLY",                      "Axoneme assembly (GO:BP)",            "Cilia",
  "GOBP_CILIUM_MOVEMENT",                       "Cilium movement (GO:BP)",             "Cilia",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION", "Translation initiation (Reactome)",   "Translation/MYC",
  "GOBP_RIBOSOME_BIOGENESIS",                   "Ribosome biogenesis (GO:BP)",         "Translation/MYC",
  "HALLMARK_MYC_TARGETS_V1",                    "MYC targets V1 (Hallmark)",           "Translation/MYC",
  "HALLMARK_MTORC1_SIGNALING",                  "mTORC1 signalling (Hallmark)",        "Translation/MYC",
  "HALLMARK_E2F_TARGETS",                       "E2F targets (Hallmark)",              "Cell cycle",
  "HALLMARK_G2M_CHECKPOINT",                    "G2/M checkpoint (Hallmark)",          "Cell cycle",
  "REACTOME_DNA_REPLICATION",                   "DNA replication (Reactome)",          "Cell cycle",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",         "Interferon-α response (Hallmark)", "IFN†",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",         "Interferon-γ response (Hallmark)", "IFN†") |>
  filter(set %in% names(all_sets))

long <- map_dfr(c(neuron, npc), load_de)
by_symbol <- function(df, value) {
  df |> filter(!is.na(Gene_Symbol)) |> group_by(Gene_Symbol) |>
    slice_max(abs({{ value }}), n = 1, with_ties = FALSE) |> ungroup() |>
    select(Gene_Symbol, v = {{ value }}) |> deframe() |> sort(decreasing = TRUE)
}
stouffer <- function(an) long |> filter(dataset %in% an) |>
  mutate(z = sign(log2FoldChange) * qnorm(pmax(pvalue, 1e-300) / 2, lower.tail = FALSE),
         w = sqrt(donors[dataset])) |>
  group_by(Gene_ID, Gene_Symbol) |> filter(n() == length(an)) |>
  summarise(z = sum(w * z) / sqrt(sum(w^2)), .groups = "drop") |> by_symbol(z)
rankings <- c(set_names(map(c(neuron, npc), ~ by_symbol(load_de(.x), stat)), c(neuron, npc)),
              list(pooled_neuron = stouffer(neuron), pooled_npc = stouffer(npc)))

dot <- imap_dfr(rankings, function(st, nm) {
  sets <- all_sets[themes$set]
  fg <- fgsea(sets, st, minSize = 10, maxSize = 500, eps = 0) |> as_tibble() |>
    transmute(set = pathway, NES, pval)
  cam <- cameraPR(st, ids2indices(sets, names(st)), inter.gene.cor = 0.05, sort = FALSE) |>
    as_tibble(rownames = "set") |>
    transmute(set, camera_dir = Direction, camera_p = PValue,
              camera_fdr = p.adjust(PValue, "BH"))
  left_join(fg, cam, by = "set") |> mutate(ranking = nm)
}) |>
  left_join(themes, by = "set") |>
  mutate(stage = if_else(ranking %in% c(neuron, "pooled_neuron"), "Neurons", "NPC"),
         stage = factor(stage, levels = c("Neurons", "NPC")),
         column = case_when(str_starts(ranking, "pooled") ~ "Pooled",
                            TRUE ~ str_extract(ranking, "^GSE\\d+")),
         column = factor(column, levels = c("GSE67528", "GSE124308", "GSE142670", "Pooled")),
         label = factor(label, levels = rev(themes$label)),
         theme = factor(theme, levels = unique(themes$theme)),
         camera_sig = camera_fdr < 0.05)
write_csv(dot, file.path(tab_dir, "Fig3_enrichment_fgsea_cameraPR.csv"))
message("cameraPR FDR < 0.05 (within each column, 12 sets):")
print(dot |> filter(camera_sig) |> select(ranking, label, camera_dir, camera_fdr, NES),
      n = Inf)

p3 <- ggplot(dot, aes(column, label)) +
  # background tile on every cell (transparent except Pooled) so the x order follows the
  # factor levels; a Pooled-only first layer would push Pooled to the front
  geom_tile(aes(alpha = column == "Pooled"), fill = "grey90", width = 1, height = 1) +
  scale_alpha_manual(values = c(`TRUE` = 0.6, `FALSE` = 0), guide = "none") +
  geom_point(aes(fill = NES, size = pmin(-log10(pval), 10), colour = camera_sig),
             shape = 21, stroke = 0.9) +
  facet_grid(theme ~ stage, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b", midpoint = 0,
                       name = "NES") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "transparent"),
                      labels = c(`TRUE` = "cameraPR\nFDR < 0.05", `FALSE` = "not significant"),
                      name = NULL, drop = FALSE) +
  scale_size_continuous(name = expression(-log[10]~italic(p)~"(fgsea)"), range = c(1, 6)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text.y = element_text(angle = 270),
        strip.background = element_rect(fill = "grey90"),
        panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4, fill = "grey80", stroke = 1),
                               order = 3))
save_fig(p3, "Fig3_enrichment_dotplot", 7.5, 5.8)

# ================================================================================
# Fig S1 — PCA panel
# ================================================================================
message("\n== Fig S1: PCA panel ==")
score_genes <- list(prolif = unique(c(all_sets$HALLMARK_E2F_TARGETS,
                                      all_sets$HALLMARK_MYC_TARGETS_V1)),
                    motile = motile)
pca_one <- function(label, gse, stage) {
  x <- read_csv(paste0("data/raw_counts_primary/", gse, "_raw_counts.csv"),
                show_col_types = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  meta <- read_csv(paste0("data/metadata/", gse, "_metadata_full.csv"),
                   show_col_types = FALSE) |> mutate(donor = as.character(donor))
  if (!is.na(stage)) meta <- filter(meta, stage == !!stage)
  meta <- filter(meta, sample %in% colnames(m))
  dge <- DGEList(m[, meta$sample])
  dge <- calcNormFactors(dge[filterByExpr(dge, group = meta$condition), , keep.lib.sizes = FALSE])
  lc  <- cpm(dge, log = TRUE, prior.count = 2)
  top <- head(order(apply(lc, 1, var), decreasing = TRUE), 500)
  pc  <- prcomp(t(lc[top, ]), scale. = FALSE)
  ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)

  # per-sample gene-set scores (mean z of expressed members)
  sym <- sym_of[rownames(lc)]
  zrow <- t(scale(t(lc)))
  sc <- map_dfc(score_genes, function(g) {
    idx <- which(sym %in% g & rowMeans(lc) > 1)
    colMeans(zrow[idx, , drop = FALSE], na.rm = TRUE)
  })
  df <- as_tibble(pc$x[, 1:2]) |> bind_cols(meta, sc) |>
    mutate(condition = factor(condition, levels = c("Control", "ASD")),
           panel = pretty[[label]],
           shape_var = if ("timepoint" %in% names(meta)) timepoint else "all")
  if (gse == "GSE124308" && "ipsc_source" %in% names(df)) {
    message("GSE124308 donors by iPSC source, with mean PC2:")
    print(df |> group_by(donor, condition, ipsc_source) |>
            summarise(mean_PC2 = round(mean(PC2), 1), .groups = "drop") |>
            arrange(desc(mean_PC2)), n = Inf)
  }
  cors <- tibble(dataset = label, PC1_var = ve[1], PC2_var = ve[2],
                 r_PC1_prolif = cor(df$PC1, df$prolif), r_PC2_prolif = cor(df$PC2, df$prolif),
                 r_PC1_motile = cor(df$PC1, df$motile), r_PC2_motile = cor(df$PC2, df$motile))
  p <- ggplot(df, aes(PC1, PC2, colour = condition)) +
    geom_line(aes(group = interaction(donor, shape_var)), linewidth = 0.3, alpha = 0.5,
              show.legend = FALSE) +
    geom_point(aes(shape = shape_var), size = 2.2) +
    geom_text_repel(aes(label = donor), size = 2.2, show.legend = FALSE, max.overlaps = 30,
                    segment.size = 0.2) +
    scale_colour_manual(values = col_cond, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c(all = 16, d35 = 17, d135 = 16), name = NULL,
                       breaks = c("d35", "d135")) +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
         subtitle = pretty[[label]]) +
    theme_bw(base_size = 9) +
    theme(plot.subtitle = element_text(face = "bold"), panel.grid.minor = element_blank())
  list(plot = p, cors = cors)
}
pcas <- pmap(tribble(~label, ~gse, ~stage,
                     "GSE67528_Neuron", "GSE67528", "Neuron",
                     "GSE124308_both", "GSE124308", NA,
                     "GSE142670_Neuron", "GSE142670", "Neuron",
                     "GSE67528_NPC", "GSE67528", "NPC",
                     "GSE142670_NPC", "GSE142670", "NPC"), pca_one)
pc_cor <- map_dfr(pcas, "cors")
message("Correlation of PC1/PC2 with proliferation (E2F + MYC V1) and motile-cilia scores:")
print(pc_cor |> mutate(across(starts_with("r_"), \(x) round(x, 2))), width = Inf)
write_csv(pc_cor, file.path(tab_dir, "FigS1_PC_score_correlations.csv"))

pS1 <- wrap_plots(c(map(pcas, "plot"), list(guide_area())), ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "A")
save_fig(pS1, "FigS1_PCA_panel", 11, 7.5)

# ================================================================================
# Fig S2 — heatmap of top-40 pooled genes
# ================================================================================
message("\n== Fig S2: heatmap ==")
pooled <- long |> filter(dataset %in% neuron) |>
  mutate(z = sign(log2FoldChange) * qnorm(pmax(pvalue, 1e-300) / 2, lower.tail = FALSE),
         w = sqrt(donors[dataset])) |>
  group_by(Gene_ID) |> filter(n() == 3) |>
  summarise(pooled_z = sum(w * z) / sqrt(sum(w^2)),
            same_sign = abs(sum(sign(log2FoldChange))) == 3,
            stars_124308 = sum(pvalue < 0.05 & dataset == "GSE124308_both"),
            .groups = "drop") |>
  filter(same_sign) |> arrange(desc(abs(pooled_z))) |> head(40)

hm <- expand_grid(Gene_ID = pooled$Gene_ID, dataset = c(neuron, npc)) |>
  left_join(long |> select(Gene_ID, dataset, log2FoldChange, pvalue),
            by = c("Gene_ID", "dataset")) |>
  left_join(pooled |> select(Gene_ID, pooled_z), by = "Gene_ID") |>
  mutate(symbol = coalesce(unname(sym_of[Gene_ID]), Gene_ID),
         star = if_else(!is.na(pvalue) & pvalue < 0.05, "*", ""),
         fill = pmax(pmin(log2FoldChange, 2), -2),
         dataset = factor(pretty[dataset], levels = pretty[c(neuron, npc)]),
         symbol = fct_reorder(symbol, pooled_z),
         cilia = symbol %in% cilia_any)
write_csv(hm |> select(symbol, Gene_ID, dataset, log2FoldChange, pvalue, pooled_z, cilia),
          file.path(tab_dir, "FigS2_heatmap_top40.csv"))
stars <- hm |> filter(star == "*") |> count(dataset, name = "genes_p_lt_0.05")
message("Genes with p < 0.05 per column (of 40):"); print(stars)
message("Cilia-annotated genes in the top 40: ", sum(distinct(hm, symbol, cilia)$cilia))

pS2 <- ggplot(hm, aes(dataset, symbol, fill = fill)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = star), size = 3, vjust = 0.75) +
  geom_point(data = distinct(hm, symbol, cilia) |> filter(cilia),
             aes(x = 0.35, y = symbol), inherit.aes = FALSE,
             shape = 15, size = 2, colour = "#d95f02") +
  geom_vline(xintercept = 3.5, linewidth = 0.6) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2),
                       labels = c("≤ -2", "-1", "0", "1", "≥ 2"),
                       na.value = "grey75", name = "log2FC\n(ASD vs control)") +
  scale_x_discrete(expand = expansion(add = c(0.9, 0.5))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
save_fig(pS2, "FigS2_heatmap_top40", 6.5, 8.5)

# ================================================================================
# Fig S3 — UpSet with expected counts for every intersection
# ================================================================================
message("\n== Fig S3: UpSet ==")
de_n <- map(set_names(neuron), load_de)
shared <- reduce(map(de_n, "Gene_ID"), intersect)
top_tbl <- map(de_n, ~ .x |> filter(Gene_ID %in% shared) |> arrange(pvalue) |> head(500))
key_of <- function(d) paste0(d$Gene_ID, if_else(d$log2FoldChange > 0, "_up", "_down"))
top_signed <- map(top_tbl, key_of)
up_share <- map_dbl(top_tbl, ~ mean(.x$log2FoldChange > 0))

combos <- list(c(1, 2), c(1, 3), c(2, 3), c(1, 2, 3))
count_combo <- function(sets) map_int(combos, ~ length(reduce(sets[.x], intersect)))
obs <- count_combo(top_signed)
perm <- replicate(2000, {
  r <- map2(top_signed, up_share, function(s, u) {
    g <- sample(shared, 500)
    paste0(g, if_else(runif(500) < u, "_up", "_down"))
  })
  count_combo(r)
})
up_tbl <- tibble(
  intersection = map_chr(combos, ~ paste(pretty[neuron][.x], collapse = " ∩ ")),
  observed = obs, expected = rowMeans(perm),
  ratio = obs / rowMeans(perm),
  perm_p_greater = map_dbl(seq_along(combos), ~ (sum(perm[.x, ] >= obs[.x]) + 1) /
                             (ncol(perm) + 1)),
  universe = length(shared))
message("Same-direction overlaps of the top 500 genes (universe = ", length(shared), "):")
print(up_tbl, width = Inf)
write_csv(up_tbl, file.path(tab_dir, "FigS3_upset_observed_expected.csv"))

names(top_signed) <- pretty[neuron]
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf_dev(file.path(fig_dir, "FigS3_upset_top500.pdf"), width = 7, height = 4.5)
  else png(file.path(fig_dir, "FigS3_upset_top500.png"), width = 7, height = 4.5,
           units = "in", res = 300)
  print(upset(fromList(top_signed), order.by = "freq", nsets = 3,
              mainbar.y.label = "Genes in top 500, same direction",
              sets.x.label = "Top-500 genes per dataset"))
  grid::grid.text(paste(sprintf("%s: %d obs / %.1f exp", c("A∩B", "A∩C", "B∩C", "A∩B∩C"),
                                obs, rowMeans(perm)), collapse = "\n"),
                  x = 0.97, y = 0.93, just = c("right", "top"),
                  gp = grid::gpar(fontsize = 7, col = "grey25"))
  dev.off()
}
message("  saved FigS3_upset_top500 (.pdf, .png); A/B/C = ", paste(pretty[neuron], collapse = " / "))

# ================================================================================
# Fig S4 — STRING network (exploratory)
# ================================================================================
message("\n== Fig S4: STRING network ==")
le_file <- "results/tables/cilia/cilia_leading_edge_pooled_neuron.csv"
if (file.exists(le_file)) {
  le <- read_csv(le_file, show_col_types = FALSE) |> arrange(pooled_z) |> head(60) |>
    mutate(class = recode(class, primary = "IFT / transition zone (all cilia)",
                          motile = "Motile-cilia specific",
                          `other/shared` = "Other cilia-annotated"))
  resp <- tryCatch(POST("https://string-db.org/api/tsv/network",
                        body = list(identifiers = paste(le$Gene_Symbol, collapse = "\r"),
                                    species = 9606, required_score = 700,
                                    caller_identity = "asd_ipsc_meta"),
                        encode = "form", timeout(120)),
                   error = function(e) NULL)
  if (!is.null(resp) && !http_error(resp)) {
    net <- read_tsv(content(resp, "text", encoding = "UTF-8"), show_col_types = FALSE)
    edges <- net |> distinct(preferredName_A, preferredName_B, score) |>
      filter(preferredName_A %in% le$Gene_Symbol, preferredName_B %in% le$Gene_Symbol)
    v <- le |> filter(Gene_Symbol %in% c(edges$preferredName_A, edges$preferredName_B))
    g <- graph_from_data_frame(edges, directed = FALSE, vertices = v)
    message("Network: ", vcount(g), " of ", nrow(le), " input genes have >= 1 edge; ",
            ecount(g), " interactions (STRING score >= 0.7)")
    write_csv(edges, file.path(tab_dir, "FigS4_string_edges.csv"))
    pS4 <- ggraph(g, layout = "stress") +
      geom_edge_link(aes(alpha = score), colour = "grey55", show.legend = FALSE) +
      geom_node_point(aes(colour = class, size = abs(pooled_z))) +
      geom_node_text(aes(label = name), size = 2.4, repel = TRUE, max.overlaps = Inf,
                     segment.size = 0.2) +
      scale_colour_manual(values = c(`Motile-cilia specific` = "#d95f02",
                                     `IFT / transition zone (all cilia)` = "#1b9e77",
                                     `Other cilia-annotated` = "grey50"), name = NULL) +
      scale_size_continuous(range = c(1.5, 5), name = "|pooled z|") +
      theme_void(base_size = 9) + theme(legend.position = "right")
    save_fig(pS4, "FigS4_string_network_exploratory", 8, 6.5)
  } else message("STRING API not reachable; network skipped.")
} else message("Run 05_cilia_check.R first; network skipped.")

# ================================================================================
# Fig S5 — volcano of the pooled neuron random-effects meta-analysis
# ================================================================================
message("\n== Fig S5: meta-analysis volcano ==")
rem_file <- "results/tables/meta/neuron_primary_REM.csv"
if (file.exists(rem_file)) {
  rem <- read_csv(rem_file, show_col_types = FALSE) |>
    filter(is.finite(est), is.finite(p)) |>
    mutate(group = case_when(Gene_Symbol %in% motile ~ "Motile-cilia gene",
                             TRUE ~ "Other"),
           group = factor(group, levels = c("Other", "Motile-cilia gene")))
  thr <- if (any(rem$fdr < 0.05)) max(rem$p[rem$fdr < 0.05]) else NA
  lab <- rem |> filter(group == "Motile-cilia gene") |> slice_min(p, n = 12)
  pS5 <- ggplot(rem, aes(est, -log10(p))) +
    geom_point(data = ~ filter(.x, group == "Other"), colour = "grey70", size = 0.5,
               alpha = 0.5) +
    geom_point(data = ~ filter(.x, group != "Other"), colour = "#d95f02", size = 1.3) +
    geom_text_repel(data = lab, aes(label = Gene_Symbol), size = 2.3, colour = "#d95f02",
                    max.overlaps = Inf, segment.size = 0.2) +
    { if (!is.na(thr)) geom_hline(yintercept = -log10(thr), linetype = 2) } +
    labs(x = "Pooled log2FC (random effects, ASD vs control)",
         y = expression(-log[10]~italic(p))) +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
  save_fig(pS5, "FigS5_volcano_meta_neurons", 5.5, 4.5)
  message("Genes tested: ", nrow(rem), "; smallest p: ", signif(min(rem$p), 2),
          "; FDR < 0.05: ", sum(rem$fdr < 0.05),
          "; motile-cilia genes with est < 0: ",
          sum(rem$est[rem$group != "Other"] < 0), "/", sum(rem$group != "Other"))
} else message("Run 02_meta_analysis.R first; volcano skipped.")

message("\nDone. Figures in ", fig_dir)
