# ==============================================================================
# 00f_get_author_counts.R — DOWNLOAD THE AUTHORS' PER-SAMPLE COUNT FILES
# ------------------------------------------------------------------------------
# GSE142670 and GSE61476 have no series-level count table; each sample (GSM)
# carries its own count file. This script:
#   1. downloads the count file(s) of every sample
#   2. shows the first lines of one file per dataset (gene ID type, columns)
#   3. if files are simple two-column tables (gene, count), assembles one
#      count matrix per dataset in the pilot format (GeneID + one column/GSM)
#
# These are the AUTHORS' counts (different annotations): used for a first look
# and for the protocol's "author counts" sensitivity analysis. The main analysis
# uses uniform Salmon re-quantification from FASTQ.
#
# Usage: pixi run Rscript scripts/00f_get_author_counts.R GSE142670 GSE61476
# Output: data/raw/author_suppl/<GSE>/<GSM>/...   (downloaded files)
#         data/raw_counts_author/<GSE>_raw_counts.csv  (if assembled)
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

options(timeout = 600)
geo_ids <- commandArgs(trailingOnly = TRUE)
if (length(geo_ids) == 0) geo_ids <- c("GSE142670", "GSE61476")

dir.create("data/raw_counts_author", showWarnings = FALSE, recursive = TRUE)

# Keep count files only (skip e.g. RPKM files)
count_regex <- "count|htseq|\\.out"

read_two_col <- function(path) {
  x <- read_tsv(path, col_names = FALSE, show_col_types = FALSE,
                comment = "#", progress = FALSE)
  if (ncol(x) != 2) return(NULL)
  x |> rename(GeneID = X1, count = X2) |>
    mutate(GeneID = str_remove_all(GeneID, "\"")) |>
    filter(!str_starts(GeneID, "__"))          # HTSeq summary rows
}

# Strip Ensembl version (".10") and PAR_Y tags, then sum any duplicate rows
clean_ids <- function(x) {
  x |>
    mutate(GeneID = str_remove(GeneID, "\\.\\d+(_PAR_Y)?$")) |>
    group_by(GeneID) |>
    summarise(count = sum(count), .groups = "drop")
}

for (geo_id in geo_ids) {
  message("\n==================== ", geo_id, " ====================")
  meta <- read_csv(file.path("data/metadata/geo_full", paste0(geo_id, "_geo_full.csv")),
                   show_col_types = FALSE)
  base <- file.path("data/raw/author_suppl", geo_id)
  dir.create(base, showWarnings = FALSE, recursive = TRUE)

  # 1. Download each sample's count file (skips samples already downloaded)
  files <- map(meta$geo_accession, function(gsm) {
    have <- list.files(file.path(base, gsm), full.names = TRUE)
    if (length(have) > 0) return(have)
    res <- tryCatch(
      getGEOSuppFiles(gsm, makeDirectory = TRUE, baseDir = base,
                      filter_regex = count_regex),
      error = function(e) { message("  ", gsm, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) character(0) else rownames(res)
  })
  names(files) <- meta$geo_accession
  n_files <- lengths(files)
  message("Samples with a count file: ", sum(n_files > 0), " / ", length(files),
          " | samples with >1 file: ", sum(n_files > 1))

  # 2. Preview one file
  first <- files[[which(n_files > 0)[1]]][1]
  message("\nFirst lines of ", basename(first), ":")
  con <- gzfile(first); writeLines(head(readLines(con, n = 8), 8)); close(con)

  # 3. Samples without a count file: which stage, and do they have other files?
  missing <- names(files)[n_files == 0]
  if (length(missing) > 0) {
    message("\nSamples with no count file: ", length(missing))
    print(meta |> filter(geo_accession %in% missing) |>
            select(geo_accession, title, source_name_ch1), n = Inf)
    message("Any other supplementary files for these samples (not downloaded):")
    walk(missing, function(gsm) {
      other <- tryCatch(getGEOSuppFiles(gsm, fetch_files = FALSE),
                        error = function(e) NULL)
      message("  ", gsm, ": ",
              if (is.null(other) || NROW(other) == 0) "none" else
                paste(other$fname, collapse = ", "))
    })
    write_csv(tibble(sample = missing),
              file.path("data/raw_counts_author", paste0(geo_id, "_missing_counts.csv")))
  }

  # 4. Assemble from samples that have exactly one two-column count file
  ok <- names(files)[n_files == 1]
  tabs <- map(files[ok], ~ read_two_col(.x[1]))
  bad  <- map_lgl(tabs, is.null)
  if (any(bad)) {
    message("\nNot two-column tables (not assembled): ", sum(bad), " samples; see preview.")
  }
  if (sum(!bad) > 0 && mean(!bad) > 0.5) {
    # ID format per file, before cleaning (versioned vs unversioned)
    id_fmt <- imap_dfr(tabs[!bad], ~ tibble(
      sample = .y, genes = nrow(.x),
      versioned_share = round(mean(str_detect(.x$GeneID, "\\.\\d+")), 3),
      example = .x$GeneID[1]))
    message("\nGene ID format per file (genes, share with a version suffix):")
    print(count(id_fmt, genes, versioned_share, name = "files"))

    mat <- imap(tabs[!bad], ~ clean_ids(.x) |> rename(!!.y := count)) |>
      reduce(full_join, by = "GeneID")
    message("\nAssembled matrix: ", nrow(mat), " genes x ", ncol(mat) - 1, " samples")
    message("Example gene IDs: ", paste(head(mat$GeneID, 6), collapse = ", "))
    message("Missing values after joining: ", sum(is.na(mat)))
    write_csv(mat, file.path("data/raw_counts_author", paste0(geo_id, "_raw_counts.csv")))
  } else {
    message("\nCount files are not simple gene/count tables; this dataset needs ",
            "FASTQ re-quantification.")
  }
}
