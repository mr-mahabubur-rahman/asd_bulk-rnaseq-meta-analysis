# ==============================================================================
# 07_figures.R — DESCRIPTIVE FIGURES FOR THE THREE-COHORT ANALYSIS
# ------------------------------------------------------------------------------
# 1. PCA per dataset-stage (QC; shows clustering by donor)
# 2. UpSet of top-500 genes per neuron dataset, split by direction, with a
#    permutation test of the 3-way same-direction overlap (replaces a DEG Venn:
#    no gene passes FDR, so these are TOP-RANKED genes, not DEGs)
# 3. Heatmap of log2FC per dataset-stage for the 40 genes with the strongest
#    pooled neuron signal (same sign in all 3 neuron datasets)
# 4. Enrichment dot plot: selected themes x dataset (fgsea NES, -log10 p)
# 5. EXPLORATORY STRING network of cilia leading-edge genes (from 05)
#
# Output: results/figures/main/*.png, results/tables/figures/*.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(ggrepel)
  library(UpSetR)
  library(igraph)
  library(ggraph)
  library(httr)
})
set.seed(2026)
fig_dir <- "results/figures/main"
tab_dir <- "results/tables/figures"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

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

# --- 1. PCA per dataset-stage ---------------------------------------------------
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
  df  <- as_tibble(pc$x[, 1:2]) |> bind_cols(meta)
  p <- ggplot(df, aes(PC1, PC2, colour = condition)) +
    geom_point(aes(shape = if ("timepoint" %in% names(df)) timepoint else "all"), size = 2.6) +
    geom_text_repel(aes(label = donor), size = 2.6, show.legend = FALSE, max.overlaps = 30) +
    scale_colour_manual(values = c(Control = "#2166ac", ASD = "#b2182b")) +
    labs(title = pretty[[label]], subtitle = "Top 500 variable genes; labels = donor",
         x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
         shape = NULL, colour = NULL) +
    theme_bw(base_size = 10)
  ggsave(file.path(fig_dir, paste0("PCA_", label, ".png")), p, width = 6, height = 5,
         dpi = 300, bg = "white")
}
pwalk(tribble(~label, ~gse, ~stage,
              "GSE67528_Neuron", "GSE67528", "Neuron",
              "GSE124308_both", "GSE124308", NA,
              "GSE142670_Neuron", "GSE142670", "Neuron",
              "GSE67528_NPC", "GSE67528", "NPC",
              "GSE142670_NPC", "GSE142670", "NPC"), pca_one)
message("PCA plots written.")

# --- 2. UpSet of top-500 same-direction genes (neurons) -------------------------
de_n <- map(set_names(neuron), load_de)
shared <- reduce(map(de_n, "Gene_ID"), intersect)
top_signed <- map(de_n, ~ .x |> filter(Gene_ID %in% shared) |> arrange(pvalue) |>
                    head(500) |>
                    mutate(key = paste0(Gene_ID, if_else(log2FoldChange > 0, "_up", "_down"))) |>
                    pull(key))
names(top_signed) <- pretty[neuron]
obs3 <- length(reduce(top_signed, intersect))

# permutation: random 500 genes per dataset with the observed up/down share
perm <- replicate(1000, {
  r <- map(de_n, function(d) {
    up_share <- mean(d |> filter(Gene_ID %in% shared) |> arrange(pvalue) |> head(500) |>
                       pull(log2FoldChange) > 0)
    g <- sample(shared, 500)
    paste0(g, if_else(runif(500) < up_share, "_up", "_down"))
  })
  length(reduce(r, intersect))
})
upset_test <- tibble(observed_3way = obs3, expected_mean = mean(perm),
                     perm_p = (sum(perm >= obs3) + 1) / (length(perm) + 1))
message("\nTop-500 same-direction 3-way overlap (neurons):")
print(upset_test)
write_csv(upset_test, file.path(tab_dir, "upset_top500_permutation.csv"))

png(file.path(fig_dir, "UpSet_top500_neurons.png"), width = 7, height = 4.5,
    units = "in", res = 300)
print(upset(fromList(top_signed), order.by = "freq", nsets = 3,
            mainbar.y.label = "Genes in top 500, same direction",
            sets.x.label = "Top-500 genes per dataset"))
grid::grid.text(sprintf("3-way overlap: %d observed vs %.1f expected (permutation p = %.3f)",
                        obs3, mean(perm), upset_test$perm_p),
                x = 0.62, y = 0.97, gp = grid::gpar(fontsize = 8))
dev.off()

# --- 3. Heatmap of log2FC for top pooled genes -----------------------------------
long <- map_dfr(c(neuron, npc), load_de)
pooled <- long |> filter(dataset %in% neuron) |>
  mutate(z = sign(log2FoldChange) * qnorm(pmax(pvalue, 1e-300) / 2, lower.tail = FALSE),
         w = sqrt(donors[dataset])) |>
  group_by(Gene_ID) |> filter(n() == 3) |>
  summarise(pooled_z = sum(w * z) / sqrt(sum(w^2)),
            same_sign = abs(sum(sign(log2FoldChange))) == 3, .groups = "drop") |>
  filter(same_sign) |> arrange(desc(abs(pooled_z))) |> head(40)
hm <- long |> filter(Gene_ID %in% pooled$Gene_ID) |>
  mutate(symbol = coalesce(unname(sym_of[Gene_ID]), Gene_ID),
         star = if_else(pvalue < 0.05, "*", ""),
         dataset = factor(pretty[dataset], levels = pretty[c(neuron, npc)])) |>
  left_join(pooled, by = "Gene_ID") |>
  mutate(symbol = fct_reorder(symbol, pooled_z))
write_csv(hm |> select(symbol, Gene_ID, dataset, log2FoldChange, pvalue, pooled_z),
          file.path(tab_dir, "heatmap_top40_log2FC.csv"))
p <- ggplot(hm, aes(dataset, symbol, fill = pmax(pmin(log2FoldChange, 2), -2))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = star), size = 3, vjust = 0.75) +
  geom_vline(xintercept = 3.5, linewidth = 0.6) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       name = "log2FC\n(ASD vs\ncontrol)") +
  labs(x = NULL, y = NULL,
       title = "Top 40 genes by pooled neuron signal (same direction in all 3 cohorts)",
       subtitle = "* p < 0.05 in that dataset; none passes FDR; NPC columns shown for comparison") +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "Heatmap_top40_log2FC.png"), p, width = 7, height = 8.5,
       dpi = 300, bg = "white")
message("Heatmap written.")

# --- 4. Enrichment dot plot --------------------------------------------------------
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
themes <- tribble(
  ~set,                                         ~theme,
  "GOBP_CILIUM_ORGANIZATION",                   "Cilia",
  "GOBP_AXONEME_ASSEMBLY",                      "Cilia",
  "GOBP_CILIUM_MOVEMENT",                       "Cilia",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION", "Translation / MYC",
  "GOBP_RIBOSOME_BIOGENESIS",                   "Translation / MYC",
  "HALLMARK_MYC_TARGETS_V1",                    "Translation / MYC",
  "HALLMARK_MTORC1_SIGNALING",                  "Translation / MYC",
  "HALLMARK_E2F_TARGETS",                       "Cell cycle",
  "HALLMARK_G2M_CHECKPOINT",                    "Cell cycle",
  "REACTOME_DNA_REPLICATION",                   "Cell cycle",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",         "Interferon (pre-specified)",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",         "Interferon (pre-specified)") |>
  filter(set %in% names(all_sets))
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
  fgsea(all_sets[themes$set], st, minSize = 10, maxSize = 500, eps = 0) |> as_tibble() |>
    transmute(set = pathway, ranking = nm, NES, pval)
}) |>
  left_join(themes, by = "set") |>
  mutate(ranking = factor(pretty[ranking], levels = pretty[c(neuron, "pooled_neuron", npc, "pooled_npc")]),
         set = fct_rev(factor(str_replace_all(set, "_", " "),
                              levels = str_replace_all(themes$set, "_", " "))),
         theme = factor(theme, levels = unique(themes$theme)))
write_csv(dot, file.path(tab_dir, "enrichment_dotplot_data.csv"))
p <- ggplot(dot, aes(ranking, set, colour = NES, size = pmin(-log10(pval), 10))) +
  geom_point() +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0) +
  scale_size_continuous(name = "-log10 p\n(capped 10)", range = c(1, 6)) +
  labs(x = NULL, y = NULL, title = "Pathway themes by dataset (fgsea; ranking by moderated t / pooled z)",
       subtitle = "fgsea p-values ignore gene-gene correlation; see cameraPR results for the strict test") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0))
ggsave(file.path(fig_dir, "Enrichment_dotplot_themes.png"), p, width = 10, height = 6,
       dpi = 300, bg = "white")
message("Enrichment dot plot written.")

# --- 5. EXPLORATORY STRING network of cilia leading-edge genes -------------------
le_file <- "results/tables/cilia/cilia_leading_edge_pooled_neuron.csv"
if (file.exists(le_file)) {
  le <- read_csv(le_file, show_col_types = FALSE) |> arrange(pooled_z) |> head(60)
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
    if (nrow(edges) == 0) stop("STRING returned no interactions at score >= 0.7")
    g <- graph_from_data_frame(edges, directed = FALSE,
                               vertices = le |> filter(Gene_Symbol %in%
                                 c(edges$preferredName_A, edges$preferredName_B)))
    p <- ggraph(g, layout = "fr") +
      geom_edge_link(aes(alpha = score), colour = "grey60", show.legend = FALSE) +
      geom_node_point(aes(colour = class, size = abs(pooled_z))) +
      geom_node_text(aes(label = name), size = 2.6, repel = TRUE) +
      scale_colour_manual(values = c(motile = "#d95f02", primary = "#1b9e77",
                                     `other/shared` = "grey50")) +
      labs(title = "EXPLORATORY: STRING network of cilia leading-edge genes",
           subtitle = "Top 60 by pooled neuron z (all lower in ASD); STRING score >= 0.7; not used for claims",
           size = "|pooled z|", colour = NULL) +
      theme_void(base_size = 10)
    ggsave(file.path(fig_dir, "Network_cilia_leading_edge_EXPLORATORY.png"), p,
           width = 8, height = 7, dpi = 300, bg = "white")
    message("Network written: ", vcount(g), " genes, ", ecount(g), " interactions.")
  } else message("STRING API not reachable; network skipped.")
} else message("Run 05_cilia_check.R first; network skipped.")
