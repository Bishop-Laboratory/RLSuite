# Goal: Get the study accessions for the R-loop mapping datasets
library(tidyverse)
source("misc/utils.R")

dataset <- read_csv("analyses/rmap_full_11_25.csv")
accessions <- dataset$GSM

load("misc/available_genomes.rda")

# Someone decided to make these private after we already analyzed their data...
# THe study ID for them is GSE130242
bad_acc <- paste0("GSM37343", 14:21)
accessions <- accessions[! accessions %in% bad_acc]


mapping <- get_public_run_info(accessions, available_genomes = available_genomes)
new_dataset <- mapping %>%
  select(GSM = accessions_original, study = condition) %>%
  full_join(dataset, by = c("GSM")) %>%
  mutate(study = case_when(
    GSM %in% bad_acc ~ "GSE130242"
  )) %>%
  write_csv("analyses/rmap_full_11_25_with_study.csv")




