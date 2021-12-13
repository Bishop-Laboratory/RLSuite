### Premise ###
# We are looking for the genomic regions that are dRNH-specific, S9.6-specific, and both
###############

# Load in the coverage for the consensus signal
library(tidyverse)
library(rtracklayer)

s96URL <- "https://rlbase-data.s3.amazonaws.com/misc/rlregions_S96.bw"
dRNHURL <- "https://rlbase-data.s3.amazonaws.com/misc/rlregions_dRNH.bw"
s96 <- rtracklayer::import.bw(s96URL)
dRNH <- rtracklayer::import.bw(dRNHURL)
message("DONE")

s96 <- s96 %>%
  as_tibble() %>%
  mutate(id = paste0(seqnames, "_", start, "_", end))

dRNH <- dRNH %>%
  as_tibble() %>%
  mutate(id = paste0(seqnames, "_", start, "_", end))

s96 <- s96 %>%
  select(id, S96 = score)
dRNH <- dRNH %>%
  select(id, dRNH = score)

dat <- full_join(s96, dRNH)
dat

dat <- dat %>%
  mutate(
    S96 = ifelse(is.na(S96), 0, S96),
    dRNH = ifelse(is.na(dRNH), 0, dRNH),
    S96_scale = S96 / 184,
    dRNH_scale = dRNH / 42
  )


hist(dat$dRNH_scale)
hist(dat$S96_scale)

pat <- "(.+)_(.+)_(.+)"
final <- dat %>%
  mutate(
    seqnames = gsub(id, pattern = pat, replacement = "\\1"),
    start = gsub(id, pattern = pat, replacement = "\\2"),
    end = gsub(id, pattern = pat, replacement = "\\3"),
    score = S96_scale - dRNH_scale
  ) %>%
  select(seqnames, start, end, score)










