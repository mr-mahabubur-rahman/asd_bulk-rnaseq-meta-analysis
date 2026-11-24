# ==============================================================================
# 05_cilia_check.R — IS THE REPLICATED "CILIA DOWN" SIGNAL REAL, AND WHICH CILIA?
# ------------------------------------------------------------------------------
# 04 found cilium / axoneme / microtubule-based-movement sets lower in ASD
# neurons in all 3 cohorts. fgsea p-values ignore gene-gene correlation, so:
#   1. cameraPR (limma) per dataset: competitive test that allows for
#      inter-gene correlation (0.01 = limma default; 0.05 = conservative).
#      Ribosome/MYC sets included for comparison.
#   2. Motile vs primary cilia: curated gene sets tested per dataset and pooled;
#      leading-edge genes of the GO cilia sets classified as motile / primary.
#   3. Ependymal / choroid-plexus markers (FOXJ1, TTR, ...): motile cilia in
#      neural cultures usually mean a few ependymal-like cells (composition).
#
# Input:  results/tables/DE_dream/*.csv, data/raw_counts_primary/*,
#         data/metadata/*_metadata_full.csv, data/annotation/ensembl115_genes.tsv
# Output: results/tables/cilia/*.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})
set.seed(2026)
out_dir <- "results/tables/cilia"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Curated cilia gene sets ---------------------------------------------------
motile <- c("DNAH1","DNAH2","DNAH3","DNAH5","DNAH6","DNAH7","DNAH9","DNAH10","DNAH11",
            "DNAH12","DNAI1","DNAI2","DNAI3","DNAI4","DNAL1","RSPH1","RSPH3","RSPH4A",
            "RSPH9","SPAG6","SPAG16","SPAG17","HYDIN","CFAP43","CFAP44","CCDC39","CCDC40",
            "ODAD1","ODAD2","DNAAF1","DNAAF2","DNAAF4","ZMYND10","FOXJ1","MCIDAS","GMNC",
            "CDC20B","TEKT1","TEKT2","TPPP3","CAPS","SNTN","DYNLRB2","PACRG","RIBC1",
            "RIBC2","CFAP45","CFAP52","ENKUR")
primary <- c("IFT20","IFT22","IFT27","IFT43","IFT46","IFT52","IFT57","IFT74","IFT80",
             "IFT81","IFT88","IFT122","IFT140","IFT172","WDR19","WDR35","TTC21B",
             "BBS1","BBS2","BBS4","BBS5","BBS7","BBS9","BBS10","BBS12","ARL6","ARL13B",
             "INPP5E","TMEM67","TMEM216","TMEM231","CC2D2A","RPGRIP1L","NPHP1","NPHP4",
             "CEP290","AHI1","TCTN1","TCTN2","MKS1","B9D1","B9D2","KIF3A","KIF3B",
             "KIFAP3","DYNC2H1","DYNC2LI1","TULP3","GPR161","SMO","ADCY3","OFD1","CEP164")
ependymal_markers <- c("FOXJ1", "TTR", "AQP1", "CLIC6", "KCNJ13", "S100B")

# --- MSigDB sets to re-test ------------------------------------------------------
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
cilia_go <- c("GOBP_CILIUM_ORGANIZATION", "GOBP_CILIUM_MOVEMENT", "GOBP_AXONEME_ASSEMBLY",
              "GOBP_MICROTUBULE_BASED_MOVEMENT",
              "GOBP_CILIUM_OR_FLAGELLUM_DEPENDENT_CELL_MOTILITY")
compare <- c("REACTOME_EUKARYOTIC_TRANSLATION_INITIATION", "GOBP_RIBOSOME_BIOGENESIS",
             "HALLMARK_MYC_TARGETS_V1")
test_sets <- c(all_sets[intersect(c(cilia_go, compare), names(all_sets))],
               list(CURATED_MOTILE_CILIA = motile, CURATED_PRIMARY_CILIA = primary))

# --- Rankings --------------------------------------------------------------------
donors <- read_csv("results/tables/DE_dream/summary.csv", show_col_types = FALSE) |>
  select(analysis, donors) |> deframe()
load_de <- function(a) read_csv(file.path("results/tables/DE_dream", paste0(a, ".csv")),
                                show_col_types = FALSE) |>
  filter(!is.na(Gene_Symbol), is.finite(stat), is.finite(pvalue)) |> mutate(dataset = a)
by_symbol <- function(df, value) {
  df |> group_by(Gene_Symbol) |> slice_max(abs({{ value }}), n = 1, with_ties = FALSE) |>
    ungroup() |> select(Gene_Symbol, v = {{ value }}) |> deframe() |> sort(decreasing = TRUE)
}
stouffer_rank <- function(analyses) {
  map_dfr(analyses, load_de) |>
    mutate(z = sign(log2FoldChange) * qnorm(pmax(pvalue, 1e-300) / 2, lower.tail = FALSE),
           w = sqrt(donors[dataset])) |>
    group_by(Gene_ID, Gene_Symbol) |> filter(n() == length(analyses)) |>
    summarise(z = sum(w * z) / sqrt(sum(w^2)), .groups = "drop") |> by_symbol(z)
}
neuron <- c("GSE67528_Neuron", "GSE124308_both", "GSE142670_Neuron")
npc    <- c("GSE67528_NPC", "GSE142670_NPC")

# --- 1. cameraPR per dataset ------------------------------------------------------
cam <- map_dfr(c(neuron, npc), function(a) {
  st  <- by_symbol(load_de(a), stat)
  idx <- ids2indices(test_sets, names(st))
  map_dfr(c(0.01, 0.05), function(r) {
    cameraPR(st, idx, inter.gene.cor = r, sort = FALSE) |>
      as_tibble(rownames = "set") |>
      transmute(dataset = a, inter_gene_cor = r, set, n_genes = NGenes,
                direction = Direction, p = PValue)
  })
})
message("\n=== 1. cameraPR per dataset (correlation-aware; direction = in ASD) ===")
print(cam |> filter(inter_gene_cor == 0.05) |>
        mutate(cell = paste0(direction, " p=", signif(p, 2))) |>
        select(set, dataset, cell) |>
        pivot_wider(names_from = dataset, values_from = cell), width = Inf)
message("(inter-gene correlation 0.05 shown; 0.01 results are in the CSV)")
write_csv(cam, file.path(out_dir, "cameraPR_per_dataset.csv"))

# --- 2. Motile vs primary: fgsea per dataset and pooled ----------------------------
curated <- list(CURATED_MOTILE_CILIA = motile, CURATED_PRIMARY_CILIA = primary)
rankings <- c(list(pooled_neuron = stouffer_rank(neuron), pooled_npc = stouffer_rank(npc)),
              set_names(map(c(neuron, npc), ~ by_symbol(load_de(.x), stat)), c(neuron, npc)))
mp <- imap_dfr(rankings, function(st, nm) {
  fgsea(curated, st, minSize = 5, maxSize = 500, eps = 0) |> as_tibble() |>
    transmute(ranking = nm, set = pathway, size, NES, p = pval)
})
message("\n=== 2a. Motile vs primary cilia (fgsea) ===")
print(mp |> mutate(across(c(NES, p), \(x) signif(x, 3))), n = Inf, width = Inf)
write_csv(mp, file.path(out_dir, "motile_vs_primary_fgsea.csv"))

# Leading edge of the GO cilia sets in the pooled neuron ranking
fg <- fgsea(all_sets[intersect(cilia_go, names(all_sets))], rankings$pooled_neuron,
            minSize = 10, maxSize = 500, eps = 0)
le <- unique(unlist(fg$leadingEdge))
le_tbl <- tibble(Gene_Symbol = le, pooled_z = rankings$pooled_neuron[le]) |>
  mutate(class = case_when(Gene_Symbol %in% motile ~ "motile",
                           Gene_Symbol %in% primary ~ "primary",
                           TRUE ~ "other/shared")) |>
  arrange(pooled_z)
message("\n=== 2b. Leading-edge genes of the GO cilia sets (pooled neurons) ===")
print(count(le_tbl, class))
message("Top 30 leading-edge genes (most negative pooled z = lowest in ASD):")
print(head(le_tbl, 30), n = 30)
write_csv(le_tbl, file.path(out_dir, "cilia_leading_edge_pooled_neuron.csv"))

# --- 3. Ependymal / choroid-plexus markers ------------------------------------------
genes <- read_tsv("data/annotation/ensembl115_genes.tsv",
                  col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE)
sym_of <- setNames(genes$external_gene_name, genes$ensembl_gene_id)
stage_of <- c(GSE67528_Neuron = "Neuron", GSE124308_both = NA, GSE142670_Neuron = "Neuron",
              GSE67528_NPC = "NPC", GSE142670_NPC = "NPC")

ep <- map_dfr(names(stage_of), function(a) {
  gse <- sub("_.*$", "", a)
  x <- read_csv(paste0("data/raw_counts_primary/", gse, "_raw_counts.csv"),
                show_col_types = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  meta <- read_csv(paste0("data/metadata/", gse, "_metadata_full.csv"),
                   show_col_types = FALSE)
  if (!is.na(stage_of[[a]])) meta <- filter(meta, stage == stage_of[[a]])
  meta <- filter(meta, sample %in% colnames(m))
  lc <- cpm(calcNormFactors(DGEList(m[, meta$sample])), log = TRUE, prior.count = 2)
  rows <- which(sym_of[rownames(lc)] %in% ependymal_markers)
  as_tibble(t(lc[rows, , drop = FALSE]), rownames = "sample") |>
    rename_with(~ unname(sym_of[.x]), -sample) |>
    left_join(meta |> mutate(donor = as.character(donor)) |>
                select(sample, donor, condition), by = "sample") |>
    pivot_longer(-c(sample, donor, condition), names_to = "gene", values_to = "logCPM") |>
    mutate(analysis = a)
})
ep_tbl <- ep |>
  group_by(analysis, gene, donor, condition) |>
  summarise(logCPM = mean(logCPM), .groups = "drop") |>
  group_by(analysis, gene) |>
  summarise(mean_logCPM = mean(logCPM),
            asd_minus_ctrl = mean(logCPM[condition == "ASD"]) -
                             mean(logCPM[condition == "Control"]),
            p_donor_t = tryCatch(t.test(logCPM ~ condition)$p.value,
                                 error = function(e) NA_real_),
            .groups = "drop")
message("\n=== 3. Ependymal / choroid-plexus markers (donor means; logCPM) ===")
message("mean_logCPM below ~0 = barely expressed (few or no such cells)")
print(ep_tbl |> mutate(across(where(is.double), \(x) signif(x, 2))), n = Inf, width = Inf)
write_csv(ep_tbl, file.path(out_dir, "ependymal_markers.csv"))
