library(tidyverse)

rmap <- read_csv("analyses/scored_rmap_full_11_17_2020.csv")
rmap %>%
  filter(mode_group %in% "DRIP",
         genome == "hg38",
         Condition == "S9.6",
         MACS2__total_peaks > 1500 & xgbverdict == 1) %>%
  mutate(peak = paste0("analyses/RLoop_Table_Analysis/data/peaks/", sample_name, "_hg38.unstranded.broadPeak")) %>%
  filter(file.exists(peak)) %>%
  pull(peak) -> peaks


system("rm -rf analyses/RLoop_Table_Analysis/data/peaksFinal"); dir.create("analyses/RLoop_Table_Analysis/data/peaksFinal", showWarnings = FALSE)

sapply(peaks, file.copy, to="analyses/RLoop_Table_Analysis/data/peaksFinal")

