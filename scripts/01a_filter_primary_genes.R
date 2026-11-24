# ==============================================================================
# 01a_filter_primary_genes.R — REMOVE ALT-HAPLOTYPE / PATCH GENE IDs  (v2)
# ------------------------------------------------------------------------------
# Why: the count matrices contain Ensembl genes annotated to alternate
# haplotype and patch contigs (e.g. CHR_HSCHR6_MHC_*), alongside the
# primary-assembly copy of the same gene. Reads are split between copies
# according to each donor's haplotype, producing spurious "DEGs"
# (DDR1, HLA-A, RING1, TRIM26, ...). This script keeps only genes on the
# primary chromosomes (1-22, X, Y, MT).
#
# v2: uses a gene table extracted from the Ensembl release-115 GTF
# (data/annotation/ensembl115_genes.tsv) instead of biomaRt, which failed
# with an SSL handshake error. A fixed Ensembl release is also better for
# reproducibility (state "Ensembl release 115" in the methods).
#
# Input:  data/annotation/ensembl115_genes.tsv  (made by the bash step)
# Output: data/annotation/gene_classification.csv
#         data/raw_counts_primary/<GSE>_raw_counts.csv
# ==============================================================================

library(tidyverse)
library(rio)

dir.create("data/raw_counts_primary", showWarnings = FALSE, recursive = TRUE)

# Datasets: from the command line, e.g.
#   pixi run Rscript scripts/01a_filter_primary_genes.R GSE142670
# default = the two pilot datasets
geo_ids <- commandArgs(trailingOnly = TRUE)
if (length(geo_ids) == 0) geo_ids <- c("GSE67528", "GSE124308")
primary_chr <- c(as.character(1:22), "X", "Y", "MT")

# Step 1: load the Ensembl gene table and classify each gene by location
anno <- read_tsv("data/annotation/ensembl115_genes.tsv",
                 col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE) |>
  # In the Ensembl GTF, unplaced/unlocalised scaffolds are named by accession
  # (KI270xxx.1, GL000xxx.1); everything else off the primary chromosomes is an
  # ALT haplotype or patch contig (e.g. HSCHR6_MHC_COX_CTG1).
  mutate(location = case_when(
    chromosome_name %in% primary_chr                ~ "primary",
    str_detect(chromosome_name, "^(KI|GL)\\d")      ~ "unplaced_scaffold",
    TRUE                                            ~ "alt_or_patch"
  ))

message("Ensembl 115 genes by location:")
print(count(anno, location))

# Step 2: read count matrices (strip any version suffix from IDs)
# Count source: data/raw_counts/ (pilot) or data/raw_counts_author/ (authors'
# per-sample counts assembled by 00f_get_author_counts.R)
count_path <- function(geo_id) {
  p <- c(paste0("data/raw_counts/", geo_id, "_raw_counts.csv"),
         paste0("data/raw_counts_author/", geo_id, "_raw_counts.csv"))
  p <- p[file.exists(p)]
  if (length(p) == 0) stop("No count file found for ", geo_id)
  message(geo_id, ": counts from ", p[1])
  p[1]
}
read_counts <- function(geo_id) {
  x <- import(count_path(geo_id))
  colnames(x)[1] <- "GeneID"
  x$GeneID <- sub("\\..*$", "", as.character(x$GeneID))
  n_na <- sum(!complete.cases(x))
  if (n_na > 0) {
    message("  ", n_na, " genes missing in some samples: dropped")
    x <- x[complete.cases(x), ]
  }
  x
}
counts_list <- set_names(map(geo_ids, read_counts), geo_ids)

# Step 3: classify every ID found in the count matrices
all_ids <- unique(unlist(map(counts_list, "GeneID")))
gene_class <- tibble(ensembl_gene_id = all_ids) |>
  left_join(anno, by = "ensembl_gene_id") |>
  mutate(location = replace_na(location, "not_in_ensembl115"))  # retired IDs

export(gene_class, paste0("data/annotation/gene_classification_",
                          paste(geo_ids, collapse = "_"), ".csv"))
message("\nIDs in the count matrices by location:")
print(count(gene_class, location))

primary_ids <- gene_class$ensembl_gene_id[gene_class$location == "primary"]

# Sanity check: these ALT copies drove the original hub genes and must be gone
check_ids <- c(DDR1_alt = "ENSG00000230456", HLA_A_alt = "ENSG00000206505",
               DDR1_primary = "ENSG00000204580", HLA_A_primary = "ENSG00000206503")
message("\nSanity check (alt = FALSE, primary = TRUE expected):")
print(tibble(label = names(check_ids), id = check_ids,
             kept = check_ids %in% primary_ids))

# Step 4: write filtered count matrices
summary_tbl <- map_dfr(geo_ids, function(geo_id) {
  x <- counts_list[[geo_id]]
  filtered <- x |> filter(GeneID %in% primary_ids)
  export(filtered, paste0("data/raw_counts_primary/", geo_id, "_raw_counts.csv"))
  loc <- gene_class$location[match(x$GeneID, gene_class$ensembl_gene_id)]
  tibble(dataset = geo_id, genes_in = nrow(x), kept = nrow(filtered),
         removed_alt_patch = sum(loc == "alt_or_patch"),
         removed_scaffold  = sum(loc == "unplaced_scaffold"),
         removed_retired   = sum(loc == "not_in_ensembl115"))
})
message("\nPer-dataset summary:")
print(summary_tbl)
