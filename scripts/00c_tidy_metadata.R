# ==============================================================================
# 00c_tidy_metadata.R — BUILD ANALYSIS-READY SAMPLE METADATA
# ------------------------------------------------------------------------------
# Reads the full GEO annotation saved by 00b_geo_metadata.R and writes one
# tidy table per dataset with the variables the new design needs:
#   sample, condition, donor, clone, stage (GSE67528) / timepoint (GSE124308)
# plus donor-level covariates that could confound ASD vs control.
#
# Output: data/metadata/<GSE>_metadata_full.csv
#         console: balance tables (donors per group, samples per donor/stage,
#                  covariates by group) to check for confounding
# ==============================================================================

library(tidyverse)

in_dir  <- "data/metadata/geo_full"
out_dir <- "data/metadata"

# ---------------------------- GSE67528 ----------------------------------------
# 14 donors (patient id), 33 iPSC clones, 3 stages: iPSC / NPC / neuron
g1 <- read_csv(file.path(in_dir, "GSE67528_geo_full.csv"), show_col_types = FALSE)

meta1 <- g1 |>
  transmute(
    sample    = geo_accession,
    title     = title,
    condition = if_else(`group:ch1` == "CONTROL", "Control", "ASD"),
    donor     = `patient id:ch1`,
    clone     = `iPSc clone:ch1`,
    stage     = recode(`cell type:ch1`,
                       "iPSCs"                           = "iPSC",
                       "iPSC-derived neural progenitors" = "NPC",
                       "iPSC-derived neurons"            = "Neuron"),
    age_at_biopsy = as.numeric(`age at biopsy:ch1`),
    sex       = `gender:ch1`
  ) |>
  mutate(condition = factor(condition, levels = c("Control", "ASD")),
         stage     = factor(stage, levels = c("iPSC", "NPC", "Neuron")))

write_csv(meta1, file.path(out_dir, "GSE67528_metadata_full.csv"))

message("\n==================== GSE67528 ====================")
message("Donors per group:")
print(meta1 |> distinct(donor, condition) |> count(condition))
message("Samples per stage x group:")
print(table(meta1$stage, meta1$condition))
message("Donors per stage x group:")
print(meta1 |> distinct(donor, condition, stage) |> count(stage, condition) |>
        pivot_wider(names_from = condition, values_from = n))
message("Samples per donor x stage (clones/replicates):")
print(table(meta1$donor, meta1$stage))
message("Age at biopsy by group (donor level):")
print(meta1 |> distinct(donor, condition, age_at_biopsy) |>
        group_by(condition) |>
        summarise(n = n(), mean_age = mean(age_at_biopsy),
                  min = min(age_at_biopsy), max = max(age_at_biopsy)))

# ---------------------------- GSE124308 ---------------------------------------
# 11 donors, cortical neurons at day 35 and day 135; title e.g. "ASD 110_135"
g2 <- read_csv(file.path(in_dir, "GSE124308_geo_full.csv"), show_col_types = FALSE)

meta2 <- g2 |>
  transmute(
    sample    = geo_accession,
    title     = title,
    condition = if_else(`disease state:ch1` == "healthy control", "Control", "ASD"),
    donor     = str_match(title, "^\\S+\\s+(\\d+)_")[, 2],
    timepoint = paste0("d", str_match(title, "_(\\d+)$")[, 2]),
    ipsc_source   = `iPSc source:ch1`,
    reprogramming = `reprogramming method:ch1`,
    sex       = `gender:ch1`
  ) |>
  mutate(condition = factor(condition, levels = c("Control", "ASD")),
         timepoint = factor(timepoint, levels = c("d35", "d135")))

# Cross-check: timepoint parsed from title must agree with the description field
stopifnot(all(str_detect(g2$description, str_remove(as.character(meta2$timepoint), "^d"))))

write_csv(meta2, file.path(out_dir, "GSE124308_metadata_full.csv"))

message("\n==================== GSE124308 ====================")
message("Donors per group:")
print(meta2 |> distinct(donor, condition) |> count(condition))
message("Samples per timepoint x group:")
print(table(meta2$timepoint, meta2$condition))
message("Reprogramming method by group (donor level):")
print(meta2 |> distinct(donor, condition, reprogramming) |>
        count(condition, reprogramming) |>
        pivot_wider(names_from = reprogramming, values_from = n, values_fill = 0))
message("iPSC source by group (donor level):")
print(meta2 |> distinct(donor, condition, ipsc_source) |>
        count(condition, ipsc_source) |>
        pivot_wider(names_from = ipsc_source, values_from = n, values_fill = 0))
message("Donor table:")
print(meta2 |> distinct(donor, condition, ipsc_source, reprogramming) |>
        arrange(condition, donor), n = Inf)
