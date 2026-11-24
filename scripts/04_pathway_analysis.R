# ==============================================================================
# 04_pathway_analysis.R — PATHWAY-LEVEL ANALYSIS (protocol Section 12)
# ------------------------------------------------------------------------------
# 1. fgsea on genes ranked by signed pooled Stouffer z:
#      neurons: GSE67528_Neuron + GSE124308_both + GSE142670_Neuron
#      NPC:     GSE67528_NPC + GSE142670_NPC
#    gene sets: GO Biological Process, Reactome, MSigDB Hallmark (10-500 genes)
# 2. Pathway-level replication: significant pathways re-tested in EACH dataset
#    (ranked by its own moderated t); count datasets with the same NES sign.
# 3. Composition check: same tests on the neuron-adjusted ranking.
# 4. Pre-specified hypothesis (from the pilot): interferon-alpha/gamma response
#    and MHC class I antigen presentation, primary-assembly genes only.
#
# Input:  results/tables/DE_dream/*.csv, results/tables/DE_dream_adjNeuron/*.csv
# Output: results/tables/pathways/*.csv, results/figures/pathways/*.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
})

set.seed(2026)
tab_dir <- "results/tables/pathways"
fig_dir <- "results/figures/pathways"
for (d in c(tab_dir, fig_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- Gene sets (works with old and new msigdbr) ----------------------------------
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
sets <- tryCatch(list(
  Hallmark = get_sets("H"),
  GO_BP    = get_sets("C5", "GO:BP"),
  Reactome = get_sets("C2", "CP:REACTOME")
), error = function(e) stop(
  "msigdbr could not load gene sets: ", conditionMessage(e),
  "\nIf it mentions 'msigdbdf', run once:\n",
  "  pixi run Rscript -e 'install.packages(\"msigdbdf\", repos = ",
  "\"https://igordot.r-universe.dev\")'"))
message("Gene sets loaded: ",
        paste(names(sets), lengths(sets), sep = " = ", collapse = ", "))

# --- Rankings ----------------------------------------------------------------------
donors <- read_csv("results/tables/DE_dream/summary.csv", show_col_types = FALSE) |>
  select(analysis, donors) |> deframe()

load_de <- function(dir, a) read_csv(file.path(dir, paste0(a, ".csv")),
                                     show_col_types = FALSE) |>
  filter(!is.na(Gene_Symbol), is.finite(stat), is.finite(pvalue)) |>
  mutate(dataset = a)

# one value per symbol: keep the gene with the largest |value|
by_symbol <- function(df, value) {
  df |> group_by(Gene_Symbol) |> slice_max(abs({{ value }}), n = 1, with_ties = FALSE) |>
    ungroup() |> select(Gene_Symbol, v = {{ value }}) |> deframe() |> sort(decreasing = TRUE)
}

stouffer_rank <- function(dir, analyses) {
  long <- map_dfr(analyses, ~ load_de(dir, .x))
  long |>
    mutate(z = sign(log2FoldChange) * qnorm(pmax(pvalue, 1e-300) / 2, lower.tail = FALSE),
           w = sqrt(donors[dataset])) |>
    group_by(Gene_ID, Gene_Symbol) |>
    filter(n() == length(analyses)) |>                   # genes in every dataset
    summarise(z = sum(w * z) / sqrt(sum(w^2)), .groups = "drop") |>
    by_symbol(z)
}

stages <- list(
  neuron = c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron"),
  npc    = c("GSE67528_NPC", "GSE142670_NPC")
)

run_all_sets <- function(stats) {
  imap_dfr(sets, function(s, coll) {
    fgsea(pathways = s, stats = stats, minSize = 10, maxSize = 500, eps = 0) |>
      as_tibble() |> mutate(collection = coll, .before = 1)
  }) |>
    group_by(collection) |> mutate(padj = p.adjust(pval, "BH")) |> ungroup() |>
    mutate(leadingEdge = map_chr(leadingEdge, ~ paste(head(.x, 15), collapse = ";"))) |>
    arrange(pval)
}

results <- list()
for (stage in names(stages)) {
  an <- stages[[stage]]
  message("\n==================== ", toupper(stage), " ====================")

  # 1. Pooled ranking
  stats <- stouffer_rank("results/tables/DE_dream", an)
  message("Ranked genes: ", length(stats))
  res <- run_all_sets(stats)
  write_csv(res, file.path(tab_dir, paste0(stage, "_fgsea_pooled.csv")))
  sig <- res |> filter(padj < 0.05)
  message("Pathways at FDR < 0.05: ",
          paste(names(table(sig$collection)), table(sig$collection), sep = " = ",
                collapse = ", "))

  # 2. Per-dataset NES for significant pathways (top 40 by p)
  top_paths <- sig |> slice_min(pval, n = 40, with_ties = FALSE)
  if (nrow(top_paths) == 0) {
    message("No pathway at FDR < 0.05 for ", stage, "; skipping per-dataset checks.")
    results[[stage]] <- res
    next
  }
  pick_sets <- function(names_wanted) {
    all <- unlist(unname(sets), recursive = FALSE)
    all[names(all) %in% names_wanted]
  }
  pw <- pick_sets(top_paths$pathway)
  per_ds <- map_dfr(an, function(a) {
    st <- by_symbol(load_de("results/tables/DE_dream", a), stat)
    fgsea(pathways = pw, stats = st, minSize = 10, maxSize = 500, eps = 0) |>
      as_tibble() |> transmute(pathway, dataset = a, NES, pval)
  })
  # 3. Neuron-adjusted pooled ranking, same pathways
  adj <- fgsea(pathways = pw,
               stats = stouffer_rank("results/tables/DE_dream_adjNeuron", an),
               minSize = 10, maxSize = 500, eps = 0) |>
    as_tibble() |> transmute(pathway, NES_adjNeuron = NES, pval_adjNeuron = pval)

  agree <- per_ds |>
    inner_join(top_paths |> select(pathway, NES_pooled = NES), by = "pathway") |>
    group_by(pathway) |>
    summarise(datasets_same_sign = sum(sign(NES) == sign(NES_pooled)),
              n_datasets = n(), .groups = "drop")
  consistency <- top_paths |>
    select(collection, pathway, NES_pooled = NES, padj_pooled = padj, size) |>
    left_join(agree, by = "pathway") |>
    left_join(adj, by = "pathway") |>
    arrange(padj_pooled)
  message("\nTop pathways: pooled NES, how many datasets agree, and neuron-adjusted NES")
  print(consistency |> mutate(across(where(is.double), \(x) signif(x, 3))) |> head(25),
        n = 25, width = Inf)
  write_csv(consistency, file.path(tab_dir, paste0(stage, "_top_pathways_consistency.csv")))
  write_csv(per_ds, file.path(tab_dir, paste0(stage, "_top_pathways_per_dataset.csv")))

  # Figure: per-dataset NES heatmap for the top 25 pathways
  {
    hm <- per_ds |> filter(pathway %in% head(consistency$pathway, 25)) |>
      mutate(pathway = str_trunc(str_replace_all(pathway, "_", " "), 60),
             pathway = fct_reorder(pathway, NES, .fun = mean))
    p <- ggplot(hm, aes(dataset, pathway, fill = NES)) +
      geom_tile(colour = "white") +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
      labs(x = NULL, y = NULL,
           title = paste0(stage, ": top pooled pathways, NES in each dataset")) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(file.path(fig_dir, paste0(stage, "_top_pathways_per_dataset.png")), p,
           width = 8.5, height = 7, dpi = 300, bg = "white")
  }
  results[[stage]] <- res
}

# --- 4. Pre-specified interferon / MHC class I test ------------------------------
prespec <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
             "HALLMARK_INTERFERON_GAMMA_RESPONSE",
             "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
             "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I")
all_sets <- unlist(unname(sets), recursive = FALSE)
missing_sets <- setdiff(prespec, names(all_sets))
if (length(missing_sets)) {
  message("\nNot found in this MSigDB version: ", paste(missing_sets, collapse = ", "))
  alt <- grep("MHC_CLASS_I($|_)", names(all_sets), value = TRUE)
  message("Closest available MHC class I sets: ", paste(head(alt, 5), collapse = ", "))
}
pre_sets <- all_sets[intersect(prespec, names(all_sets))]

pre_tbl <- map_dfr(names(stages), function(stage) {
  rankings <- c(list(pooled = stouffer_rank("results/tables/DE_dream", stages[[stage]])),
                set_names(map(stages[[stage]],
                              ~ by_symbol(load_de("results/tables/DE_dream", .x), stat)),
                          stages[[stage]]))
  imap_dfr(rankings, function(st, nm) {
    fgsea(pathways = pre_sets, stats = st, minSize = 10, maxSize = 500, eps = 0) |>
      as_tibble() |> transmute(stage, ranking = nm, pathway, size, NES, pval)
  })
}) |>
  group_by(stage, ranking) |> mutate(padj_within_family = p.adjust(pval, "BH")) |>
  ungroup()
message("\n=== Pre-specified interferon / MHC class I test ",
        "(BH within the family of ", length(pre_sets), " sets) ===")
print(pre_tbl |> mutate(across(where(is.double), \(x) signif(x, 3))), n = Inf, width = Inf)
write_csv(pre_tbl, file.path(tab_dir, "prespecified_ifn_mhc1.csv"))
