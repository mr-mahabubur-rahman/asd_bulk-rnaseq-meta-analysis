# ==============================================================================
# 06_cilia_foxj1_adjusted.R — DOES THE PRIMARY-CILIA SIGNAL SURVIVE ADJUSTMENT
#                             FOR FOXJ1 / MOTILE-CILIA CELL CONTENT?
# ------------------------------------------------------------------------------
# 05 showed: motile-cilia genes lower in ASD (likely FOXJ1+ ependymal-like cell
# content; FOXJ1 lower in ASD in 5/5 dataset-stages), and primary-cilia genes
# also lower. Here each model is re-fitted with a covariate for that content:
#   version "adjFOXJ1"  : + FOXJ1 expression (z-scored logCPM)
#   version "adjMotile" : + curated motile-cilia score (mean z of expressed genes)
# Primary-cilia genes are not in either covariate, so testing them after
# adjustment is not circular. Motile-set results after "adjMotile" are
# circular by design and are shown only as a check that the adjustment worked.
#
# Dataset-stages with a cilia signal in 05: 3 neuron datasets + GSE142670_NPC.
# Output: results/tables/DE_dream_adjFOXJ1/, results/tables/DE_dream_adjMotile/,
#         results/tables/cilia/cilia_adjusted_comparison.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
  library(fgsea)
  library(msigdbr)
})
set.seed(2026)
param <- MulticoreParam(workers = 4)

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
go <- get_sets("C5", "GO:BP")
test_sets <- c(list(CURATED_PRIMARY_CILIA = primary, CURATED_MOTILE_CILIA = motile),
               go[intersect(c("GOBP_CILIUM_ORGANIZATION", "GOBP_AXONEME_ASSEMBLY"), names(go))])

genes <- read_tsv("data/annotation/ensembl115_genes.tsv",
                  col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE)
sym_of <- setNames(genes$external_gene_name, genes$ensembl_gene_id)

analyses <- tribble(
  ~label,             ~gse,        ~stage,   ~base,
  "GSE67528_Neuron",  "GSE67528",  "Neuron", "~ condition + (1|donor)",
  "GSE124308_both",   "GSE124308", NA,       "~ condition + timepoint + (1|donor)",
  "GSE142670_Neuron", "GSE142670", "Neuron", "~ condition + (1|donor)",
  "GSE142670_NPC",    "GSE142670", "NPC",    "~ condition + (1|donor)"
)

load_inputs <- function(gse, stage) {
  x <- read_csv(paste0("data/raw_counts_primary/", gse, "_raw_counts.csv"),
                show_col_types = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  meta <- read_csv(paste0("data/metadata/", gse, "_metadata_full.csv"),
                   show_col_types = FALSE) |>
    mutate(condition = factor(condition, levels = c("Control", "ASD")),
           donor = factor(as.character(donor)))
  if (!is.na(stage)) meta <- filter(meta, stage == !!stage)
  meta <- filter(meta, sample %in% colnames(m)) |> as.data.frame()
  if ("timepoint" %in% names(meta)) meta$timepoint <- factor(meta$timepoint)
  rownames(meta) <- meta$sample
  list(counts = m[, meta$sample], meta = meta)
}

add_covariates <- function(counts, meta) {
  lc <- cpm(calcNormFactors(DGEList(counts)), log = TRUE, prior.count = 2)
  sym <- sym_of[rownames(lc)]
  block <- if ("timepoint" %in% names(meta)) meta$timepoint else rep("all", nrow(meta))
  z <- function(v) ave(v, block, FUN = function(x) {
    s <- sd(x); if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s })
  meta$foxj1 <- z(lc[which(sym == "FOXJ1")[1], ])
  idx <- which(sym %in% motile & rowMeans(lc) > 1)
  meta$motile <- colMeans(t(apply(lc[idx, , drop = FALSE], 1, z)))
  message("  motile score from ", length(idx), " expressed genes; ",
          "cor(FOXJ1, motile score) = ", round(cor(meta$foxj1, meta$motile), 2))
  meta
}

fit_dream <- function(counts, meta, form) {
  dge  <- DGEList(counts)
  keep <- filterByExpr(dge, group = meta$condition)
  dge  <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  vobj <- voomWithDreamWeights(dge, form, meta, BPPARAM = param)
  fit  <- eBayes(dream(vobj, form, meta, BPPARAM = param))
  tt   <- topTable(fit, coef = "conditionASD", number = Inf, sort.by = "none")
  tibble(Gene_ID = rownames(tt), log2FoldChange = tt$logFC, lfcSE = abs(tt$logFC / tt$t),
         stat = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val,
         Gene_Symbol = unname(sym_of[rownames(tt)]))
}

set_tests <- function(res, version, label) {
  st <- res |> filter(!is.na(Gene_Symbol), is.finite(stat)) |>
    group_by(Gene_Symbol) |> slice_max(abs(stat), n = 1, with_ties = FALSE) |>
    ungroup() |> select(Gene_Symbol, stat) |> deframe() |> sort(decreasing = TRUE)
  cam <- cameraPR(st, ids2indices(test_sets, names(st)), inter.gene.cor = 0.05,
                  sort = FALSE) |> as_tibble(rownames = "set")
  fg <- fgsea(test_sets, st, minSize = 5, maxSize = 500, eps = 0) |> as_tibble()
  cam |> transmute(set, camera_dir = Direction, camera_p = PValue) |>
    left_join(fg |> transmute(set = pathway, NES), by = "set") |>
    mutate(analysis = label, version = version, .before = 1)
}

out <- list()
for (i in seq_len(nrow(analyses))) {
  a <- analyses[i, ]
  message("\n==================== ", a$label, " ====================")
  inp  <- load_inputs(a$gse, a$stage)
  meta <- add_covariates(inp$counts, inp$meta)

  dt <- meta |> group_by(donor, condition) |>
    summarise(foxj1 = mean(foxj1), motile = mean(motile), .groups = "drop")
  message("  FOXJ1 (donor means) ASD - control = ",
          round(mean(dt$foxj1[dt$condition == "ASD"]) -
                mean(dt$foxj1[dt$condition == "Control"]), 2), " SD, p = ",
          signif(t.test(foxj1 ~ condition, data = dt)$p.value, 2))

  unadj <- read_csv(file.path("results/tables/DE_dream", paste0(a$label, ".csv")),
                    show_col_types = FALSE)
  out[[paste(a$label, "unadjusted")]] <- set_tests(unadj, "unadjusted", a$label)

  for (cv in c("foxj1", "motile")) {
    form <- as.formula(paste(a$base, "+", cv))
    message("  fitting ", deparse(form))
    res <- fit_dream(inp$counts, meta, form)
    dir <- paste0("results/tables/DE_dream_adj", if (cv == "foxj1") "FOXJ1" else "Motile")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    write_csv(res, file.path(dir, paste0(a$label, ".csv")))
    out[[paste(a$label, cv)]] <- set_tests(res, paste0("adj_", cv), a$label)
  }
}

cmp <- bind_rows(out)
dir.create("results/tables/cilia", showWarnings = FALSE, recursive = TRUE)
write_csv(cmp, "results/tables/cilia/cilia_adjusted_comparison.csv")

message("\n=== Primary-cilia set: unadjusted vs adjusted (cameraPR, cor 0.05; NES from fgsea) ===")
show <- function(s) cmp |> filter(set == s) |>
  mutate(cell = paste0(camera_dir, " p=", signif(camera_p, 2), " NES=", round(NES, 2))) |>
  select(analysis, version, cell) |>
  pivot_wider(names_from = version, values_from = cell)
print(show("CURATED_PRIMARY_CILIA"), width = Inf)
message("\n=== Motile-cilia set (adj_motile is circular: a check that adjustment worked) ===")
print(show("CURATED_MOTILE_CILIA"), width = Inf)
message("\n=== GO: cilium organization ===")
print(show("GOBP_CILIUM_ORGANIZATION"), width = Inf)
