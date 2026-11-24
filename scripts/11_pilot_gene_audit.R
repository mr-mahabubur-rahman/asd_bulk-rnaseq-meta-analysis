# ==============================================================================
# 11_pilot_gene_audit.R — WHAT HAPPENED TO THE 39 PILOT "META-DEGs"?
# ------------------------------------------------------------------------------
# For each gene in the withdrawn pilot list (results/tables/meta-analysis/
# filtered_meta_degs.csv; also in git tag pilot-v0):
#   1. annotation class in Ensembl 115: primary assembly, alternate haplotype /
#      patch (which region), unplaced scaffold, or not in release 115
#   2. for primary-assembly genes: the donor-aware result in every primary model
#      (FDR, p, direction vs the pilot) and in the neuron meta-analysis
# Output: results/tables/pilot_audit/pilot_39_audit.csv + a summary in the log
# Runtime: seconds
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))
out_dir <- "results/tables/pilot_audit"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pilot_file <- "results/tables/meta-analysis/filtered_meta_degs.csv"
if (!file.exists(pilot_file)) {
  # restore the pilot table from the archived tag if it is no longer on disk
  system(paste("git show pilot-v0:", pilot_file, " > /tmp/pilot_degs.csv", sep = ""))
  pilot_file <- "/tmp/pilot_degs.csv"
}
pilot <- read_csv(pilot_file, show_col_types = FALSE) |>
  transmute(Gene_ID, pilot_symbol = Gene_Symbol, pilot_log2FC = randomSummary,
            pilot_p = randomP)
message("Pilot genes: ", nrow(pilot))

genes <- read_tsv("data/annotation/ensembl115_genes.tsv",
                  col_types = cols(.default = "c"), na = "NA") |>
  distinct(ensembl_gene_id, .keep_all = TRUE)
primary_chr <- c(as.character(1:22), "X", "Y", "MT")

audit <- pilot |>
  left_join(genes |> select(Gene_ID = ensembl_gene_id, chromosome_name,
                            symbol_115 = external_gene_name), by = "Gene_ID") |>
  mutate(class = case_when(
    is.na(chromosome_name)                  ~ "not in Ensembl 115",
    chromosome_name %in% primary_chr        ~ "primary assembly",
    str_detect(chromosome_name, "^(KI|GL)\\d") ~ "unplaced scaffold",
    TRUE                                    ~ "alternate haplotype / patch"),
    mhc_region = str_detect(coalesce(chromosome_name, ""), "MHC"))

message("\n1. Annotation class of the pilot genes:")
print(count(audit, class))
message("Of the non-primary genes, on MHC alternate haplotypes: ",
        sum(audit$mhc_region), " of ", sum(audit$class != "primary assembly"))
message("Distinct gene symbols among the 39: ",
        n_distinct(coalesce(audit$symbol_115, audit$pilot_symbol)))

# --- 2. donor-aware results for the primary-assembly genes ---------------------
models <- c("GSE67528_iPSC", "GSE67528_NPC", "GSE67528_Neuron", "GSE124308_both")
de <- map_dfr(models, function(m) {
  f <- file.path("results/tables/DE_dream", paste0(m, ".csv"))
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE) |>
    select(Gene_ID, log2FoldChange, pvalue, padj) |> mutate(model = m)
})
prim <- audit |> filter(class == "primary assembly") |> select(Gene_ID, pilot_symbol, pilot_log2FC)
per_model <- prim |> inner_join(de, by = "Gene_ID") |>
  mutate(same_dir = sign(log2FoldChange) == sign(pilot_log2FC))
best <- per_model |> group_by(Gene_ID) |>
  summarise(models_tested = n(),
            min_fdr = min(padj, na.rm = TRUE),
            n_p_lt_0.05 = sum(pvalue < 0.05),
            n_same_dir = sum(same_dir), .groups = "drop")

rem_file <- "results/tables/meta/neuron_primary_REM.csv"
rem <- if (file.exists(rem_file)) read_csv(rem_file, show_col_types = FALSE) |>
  select(Gene_ID, meta_est = est, meta_p = p, meta_fdr = fdr) else tibble(Gene_ID = character())

audit <- audit |> left_join(best, by = "Gene_ID") |> left_join(rem, by = "Gene_ID")
write_csv(audit, file.path(out_dir, "pilot_39_audit.csv"))
write_csv(per_model, file.path(out_dir, "pilot_primary_genes_per_model.csv"))

message("\n2. Primary-assembly pilot genes under donor-aware models (GSE67528 stages, GSE124308):")
print(audit |> filter(class == "primary assembly") |>
        select(pilot_symbol, pilot_log2FC, models_tested, min_fdr, n_p_lt_0.05,
               n_same_dir, meta_p, meta_fdr) |>
        mutate(across(where(is.double), \(x) signif(x, 2))), n = Inf, width = Inf)
message("Primary-assembly genes with FDR < 0.05 in any donor-aware model: ",
        sum(audit$min_fdr < 0.05, na.rm = TRUE), " of ", sum(audit$class == "primary assembly"))
message("... with FDR < 0.05 in the neuron meta-analysis: ",
        sum(audit$meta_fdr < 0.05, na.rm = TRUE))
message("\nDone. Tables in ", out_dir)
