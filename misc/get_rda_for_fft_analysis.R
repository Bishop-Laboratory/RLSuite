# Script that gets the .rda files for Dan's FFT testing

library(tidyverse)

dir.create("misc/datasets_for_fft_testing", showWarnings = FALSE)

# Manually selected samples that are representative of R-loop mapping
cases <- c(
  "SRX1070676", "SRX1070678", "SRX4231163", "SRX3084730", 
  "SRX1025890", "SRX2187023", "SRX534094", "SRX255721", 
  "ERX2277510", "SRX6427717", "SRX6651643", "SRX8908657", 
  "SRX5729069", "SRX3581338", "SRX8122753", "SRX1798983",
  "SRX3084735", "SRX1025898", "SRX2733011", "SRX2642938",
  "SRX2642947", "SRX6427721", "SRX8122770", "SRX1070677",
  "SRX3084734", "SRX1025902", "SRX2733014"
)
rmapsamps <- read_csv("analyses/rmap_full_11_25.csv") %>%
  filter(genome == "hg38",
         mode %in% c("DRIP", "DRIPc", "qDRIP", "sDRIP"),
         Condition %in% c("S9.6", "RNH")) %>%
  select(SRX, sample_name, Condition, clean_name) %>%
  mutate(group = case_when(Condition == "RNH" ~ "control", 
                           SRX %in% !! cases ~ "case",
                           TRUE ~ as.character(NA))) %>%
  filter(! is.na(group)) %>%
  mutate(filename = paste0( "misc/datasets_for_fft_testing/rda_qc/", 
                            sample_name, "_hg38.QC_report.rda")) %>%
  write_csv("misc/datasets_for_fft_testing/manifest.csv")

# Keep the qual reports needed for this analysis
files <- list.files("misc/datasets_for_fft_testing/rda_qc", full.names = TRUE) 
file.remove(files[which(! files %in% rmapsamps$filename)])

# Improve space savings by using xz compression
files <- list.files("misc/datasets_for_fft_testing/rda_qc", full.names = TRUE) 
sapply(files, function(x) {
  load(x)
  save(data_list, file = x, compress = "xz")
})

# Tar it up
system("cd misc && tar c datasets_for_fft_testing/ | xz > datasets_for_fft_testing.tar.xz")

# Uploaded to AWS at https://rmapdb-data.s3.us-east-2.amazonaws.com/misc/datasets_for_fft_testing.tar.xz
