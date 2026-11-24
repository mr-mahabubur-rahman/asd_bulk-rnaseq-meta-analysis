# ==============================================================================
# 00d_inspect_geo.R — INSPECT A NEW GEO DATASET BEFORE ANALYSIS
# ------------------------------------------------------------------------------
# For each accession given on the command line:
#   1. downloads the series-matrix sample annotation and saves it
#   2. prints every annotation field (to find donor, clone, stage, diagnosis)
#   3. prints how the authors processed the data (aligner, annotation)
#   4. lists the supplementary files on GEO (are gene counts available?)
#   5. lists SRA links (raw reads for re-quantification)
#
# Usage (from the project root):
#   pixi run Rscript scripts/00d_inspect_geo.R GSE142670 GSE61476
#
# Output: data/metadata/geo_full/<GSE>_geo_full.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

options(timeout = 600)
geo_ids <- commandArgs(trailingOnly = TRUE)
if (length(geo_ids) == 0) geo_ids <- "GSE142670"

out_dir <- "data/metadata/geo_full"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

show_field <- function(vals, col) {
  if (n_distinct(vals) <= 30) {
    message("\n--- ", col, " ---")
    print(table(vals, useNA = "ifany"))
  } else {
    message("\n--- ", col, " (", n_distinct(vals), " distinct values; first 8) ---")
    print(head(unique(vals), 8))
  }
}

for (geo_id in geo_ids) {
  message("\n==================== ", geo_id, " ====================")

  # 1. Sample annotation (a series can span several platforms / sub-series)
  gse <- getGEO(geo_id, GSEMatrix = TRUE, getGPL = FALSE)
  pd  <- map_dfr(gse, ~ as_tibble(pData(phenoData(.x))))
  message("Samples: ", nrow(pd), " | platforms: ",
          paste(unique(pd$platform_id), collapse = ", "))
  write_csv(pd, file.path(out_dir, paste0(geo_id, "_geo_full.csv")))

  # 2. Sample-level fields
  fields <- pd |>
    select(title, source_name_ch1, any_of(c("description", "library_strategy",
                                            "library_source", "library_selection",
                                            "instrument_model")),
           ends_with(":ch1"))
  for (col in names(fields)) show_field(fields[[col]], col)

  # 3. How the authors processed the data (first sample; usually identical)
  proc <- pd |> select(starts_with("data_processing")) |> slice(1) |>
    unlist() |> unname()
  message("\n--- Authors' data processing (first sample) ---")
  if (length(proc)) writeLines(paste0("  * ", proc)) else message("  (none listed)")

  # 4. Supplementary files available on GEO (listed, not downloaded)
  message("\n--- Supplementary files on GEO ---")
  supp <- tryCatch(getGEOSuppFiles(geo_id, fetch_files = FALSE),
                   error = function(e) { message("  could not list: ",
                                                 conditionMessage(e)); NULL })
  if (!is.null(supp)) print(supp)

  # 5. Raw-read links (SRA / BioSample)
  rel <- pd |> select(starts_with("relation")) |> slice(1:3)
  message("\n--- Relation fields (first 3 samples) ---")
  print(rel, width = Inf)
}

message("\nDone. Share this log so we can decide donor/stage parsing and the ",
        "count source (GEO supplementary counts vs FASTQ re-quantification).")
