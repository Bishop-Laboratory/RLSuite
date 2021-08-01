library(tidyverse)
library(magrittr)

system("rm -rf analyses/RLoop_Table_Analysis/data/peaksFinal"); dir.create("analyses/RLoop_Table_Analysis/data/peaksFinal", showWarnings = FALSE)
rmap <- read_csv("analyses/scored_rmap_full_11_17_2020.csv")
rmap %<>%
  dplyr::filter(mode_group %in% "DRIP",
         genome == "hg38",
         Condition == "S9.6",
         MACS2__total_peaks > 1500 & xgbverdict == 1) %>%
  dplyr::mutate(peak = paste0("analyses/RLoop_Table_Analysis/data/peaks/", sample_name, "_hg38.unstranded.broadPeak")) %>%
  dplyr::filter(file.exists(peak)) %T>% {
    pull(., peak) %>%
    sapply(file.copy, to="analyses/RLoop_Table_Analysis/data/peaksFinal")
  }
rmap %>%
  mutate(Cell = case_when(
    Group == "Li" ~ "Breast Tissue",
    TRUE ~ Cell
  )) %>%
  dplyr::select(Sample=clean_name, Mode=mode, Tissue=Cell, `RLoops Detected`=MACS2__total_peaks)  %>%
  tableone::CreateTableOne(vars = c("RLoops Detected", "Mode", "Tissue"), data=.) %>%
  tableone::kableone()


RLFS_PEAKS <- "hg38.rlfs.bed"
RLFS_PEAKS_URI <- "s3://rmapdb-data/misc/rlfs.hg38.fixed.bed"
download.file("https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/hg38.rlfs.bed", 
              destfile = "misc/hg38.rlfs.bed")
rlfs <- rtracklayer::import("misc/hg38.rlfs.bed")
names(rlfs) <- paste0("hg38_", seq(width(rlfs)))
score(rlfs) <- 0
rtracklayer::export(rlfs, con = RLFS_PEAKS)
system(paste0("aws s3 --region us-west-2 cp ", RLFS_PEAKS, " ", RLFS_PEAKS_URI, " --acl public-read"))

