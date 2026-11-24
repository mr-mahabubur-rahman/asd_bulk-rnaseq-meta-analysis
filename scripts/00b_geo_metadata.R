# ==============================================================================
# 00b_geo_metadata.R — DOWNLOAD FULL SAMPLE ANNOTATION FROM GEO
# ------------------------------------------------------------------------------
# Why: the current metadata files only contain `sample` + `condition`.
# To model the data correctly we need donor, clone, cell stage / time point
# and sex for every sample. GEO stores these in the series-matrix
# "characteristics" fields, which this script pulls and saves.
#
# Output:
#   data/metadata/geo_full/<GSE>_geo_full.csv   (all annotation fields)
#   console: a table of values for each field, so we can see which field
#            encodes donor / stage / sex
# ==============================================================================

library(GEOquery)
library(tidyverse)

options(timeout = 300)   # GEO downloads can be slow

out_dir <- "data/metadata/geo_full"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

geo_ids <- c("GSE67528", "GSE124308")

for (geo_id in geo_ids) {

  message("\n==================== ", geo_id, " ====================")

  # Step 1: download series matrix (a list if the series has >1 platform)
  gse <- getGEO(geo_id, GSEMatrix = TRUE, getGPL = FALSE)

  # Step 2: combine phenotype data across platforms
  pd <- map_dfr(gse, ~ as_tibble(pData(phenoData(.x)), rownames = NULL))

  # Step 3: keep the informative columns
  #   - title / source_name often encode donor + clone + cell type
  #   - "<field>:ch1" columns are the parsed characteristics (e.g. "cell type:ch1")
  meta_full <- pd |>
    select(geo_accession, title, source_name_ch1,
           any_of("description"),
           ends_with(":ch1"),
           starts_with("characteristics_ch1"))

  # Step 4: check overlap with the count matrix columns
  counts_hdr <- names(read_csv(paste0("data/raw_counts/", geo_id, "_raw_counts.csv"),
                               n_max = 0, show_col_types = FALSE))[-1]
  message("Samples in GEO: ", nrow(meta_full),
          " | in count matrix: ", length(counts_hdr),
          " | overlap: ", sum(meta_full$geo_accession %in% counts_hdr))

  # Step 5: save
  write_csv(meta_full, file.path(out_dir, paste0(geo_id, "_geo_full.csv")))

  # Step 6: print the values of each field (skip near-unique fields like title)
  for (col in setdiff(names(meta_full), "geo_accession")) {
    vals <- meta_full[[col]]
    if (n_distinct(vals) <= 30) {
      message("\n--- ", col, " ---")
      print(table(vals, useNA = "ifany"))
    } else {
      message("\n--- ", col, " (", n_distinct(vals), " distinct values; first 10) ---")
      print(head(unique(vals), 10))
    }
  }
}

message("\nDone. Please share the console output (or the two CSVs) so we can ",
        "decide which fields define donor, clone, stage/time point and sex.")
