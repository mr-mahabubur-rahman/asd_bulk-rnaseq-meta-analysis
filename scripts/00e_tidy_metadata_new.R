# ==============================================================================
# 00e_tidy_metadata_new.R — ANALYSIS-READY METADATA FOR GSE142670 AND GSE61476
# ------------------------------------------------------------------------------
# Reads the GEO annotation saved by 00d_inspect_geo.R and writes one tidy table
# per dataset (same layout as the pilot: sample, condition, donor, clone, stage)
# plus family and sex where available. Prints balance tables to check parsing.
#
# GSE142670 (Griesi-Oliveira 2021): title "F2613.cl1" = donor F2613, clone 1
#   stage: NPC / mixed neuronal cells (-> Neuron) / sorted neurons (-> Neuron_sorted)
# GSE61476 (Mariani 2015): donor id "03-01" = family 03, member 01
#   description "03-03_4_11" = donor_clone_day; days 0, 11, 31
#
# Output: data/metadata/GSE142670_metadata_full.csv
#         data/metadata/GSE61476_metadata_full.csv
# ==============================================================================

library(tidyverse)

in_dir  <- "data/metadata/geo_full"
out_dir <- "data/metadata"

# ---------------------------- GSE142670 ---------------------------------------
g1 <- read_csv(file.path(in_dir, "GSE142670_geo_full.csv"), show_col_types = FALSE)

meta1 <- g1 |>
  transmute(
    sample    = geo_accession,
    title     = title,
    condition = if_else(`disease status:ch1` == "Control", "Control", "ASD"),
    donor     = str_extract(title, "^[^.]+"),           # text before the first "."
    # clone number: "cl1" / "CL1" anywhere in the title ("F2613.GFP.CL1" -> F2613.CL1)
    clone     = paste0(donor, ".CL",
                       str_extract(title, regex("(?<=cl)\\d+", ignore_case = TRUE))),
    stage     = recode(source_name_ch1,
                       "NPC"                  = "NPC",
                       "mixed neuronal cells" = "Neuron",
                       "neuron"               = "Neuron_sorted")
  ) |>
  mutate(condition = factor(condition, levels = c("Control", "ASD")),
         stage     = factor(stage, levels = c("NPC", "Neuron", "Neuron_sorted")))

message("\n==================== GSE142670 ====================")
odd <- meta1 |> filter(str_ends(clone, "NA"))
if (nrow(odd) > 0) {
  message("Titles that do not match 'DONOR.clN' (check these by hand):")
  print(odd |> select(sample, title, stage))
}
message("Donors per group:")
print(meta1 |> distinct(donor, condition) |> count(condition))
message("Every donor must have ONE diagnosis (n_dx should be 1):")
print(meta1 |> group_by(donor) |> summarise(n_dx = n_distinct(condition), .groups = "drop") |>
        count(n_dx))
message("Samples per donor x stage:")
print(table(meta1$donor, meta1$stage))
message("Donors per stage x group:")
print(meta1 |> distinct(donor, condition, stage) |> count(stage, condition) |>
        pivot_wider(names_from = condition, values_from = n))

write_csv(meta1, file.path(out_dir, "GSE142670_metadata_full.csv"))

# ---------------------------- GSE61476 ----------------------------------------
g2 <- read_csv(file.path(in_dir, "GSE61476_geo_full.csv"), show_col_types = FALSE)

meta2 <- g2 |>
  transmute(
    sample    = geo_accession,
    title     = title,
    description = description,
    condition = if_else(str_detect(`donor's disease state:ch1`, "autism"), "ASD", "Control"),
    donor     = `donor id:ch1`,
    family    = str_extract(`donor id:ch1`, "^[^-]+"),
    day       = paste0("d", `number of days in vitro in terminal differentiation conditions:ch1`),
    sex       = `donor's gender:ch1`
  ) |>
  # Two description formats:
  #   short "1123-01_1_0"        = donor _ clone _ day
  #   long  "03-01_23_2_11_005"  = donor _ run code _ clone _ day _ library no.
  mutate(
    parts  = str_split(description, "_"),
    n_part = lengths(parts),
    clone  = paste0(donor, "_cl", map2_chr(parts, n_part, ~ if (.y == 5) .x[3] else .x[2])),
    day_in_desc = map2_chr(parts, n_part, ~ if (.y == 5) .x[4] else .x[3]),
    batch  = if_else(n_part == 5, paste0("run", map_chr(parts, 2)), "short_format")
  ) |>
  select(-parts, -n_part) |>
  mutate(condition = factor(condition, levels = c("Control", "ASD")),
         day       = factor(day, levels = c("d0", "d11", "d31")))

# Cross-check: day parsed from the description must match the day field
day_from_desc <- paste0("d", meta2$day_in_desc)
message("\n==================== GSE61476 ====================")
agree <- day_from_desc == as.character(meta2$day)
message("Description format / batch by group (samples):")
print(table(meta2$batch, meta2$condition))
message("Day field agrees with description for all samples: ", all(agree))
if (!all(agree)) {
  message("Rows where they disagree (the 'day' field is used; check clone parsing):")
  print(meta2 |> filter(!agree) |> select(sample, title, description, donor, clone, day),
        n = Inf)
}
message("Clones per donor:")
print(meta2 |> distinct(donor, clone) |> count(donor, name = "clones"), n = Inf)

message("Donors per group:")
print(meta2 |> distinct(donor, condition) |> count(condition))
message("Donor table (family, diagnosis, sex):")
print(meta2 |> distinct(family, donor, condition, sex) |> arrange(family, donor), n = Inf)
message("Samples per donor x day:")
print(table(meta2$donor, meta2$day))
message("Donors per day x group:")
print(meta2 |> distinct(donor, condition, day) |> count(day, condition) |>
        pivot_wider(names_from = condition, values_from = n))

write_csv(meta2, file.path(out_dir, "GSE61476_metadata_full.csv"))
