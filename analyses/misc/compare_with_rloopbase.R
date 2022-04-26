## Compare with the supplementary table from R-loop Base
## How many of their "high-quality" samples are POS in our analysis?

library(tidyverse)
library(readxl)
library(SRAdb)

samps_us <- RLHub::rlbase_samples()

samps_them_orig <- read_excel("analyses/misc/Supplementary Tables revised.xlsx")


# Wrangle
samps_them <- samps_them_orig %>%
  mutate(run = strsplit(`SRA Identifier`, ", ")) %>%
  unnest(run)

# Get SRX
rlp <- read_tsv("analyses/misc/rlpipes_config.tsv")   # From RLBase v1.0 data pipeline
samps_them <- inner_join(samps_them, select(rlp, experiment, run), by = c("run"))

# Get verdicts
samps_us %>%
  filter(rlsample %in% samps_them$experiment) %>%
  group_by(label, prediction) %>% tally()

samps_us %>%
  filter(rlsample %in% samps_them$experiment) %>%
  filter(prediction == "NEG", label == "POS") 

final <- samps_us %>% 
  select(rlsample, label, prediction) %>%
  right_join(unique(select(samps_them, -run)), by = c("rlsample" = "experiment")) %>%
  full_join(samps_them_orig) %>%
  distinct(`SRA Identifier`, .keep_all = TRUE) %>%
  arrange(match(Experiment, samps_them_orig$Experiment)) 

final %>%
  filter(label == "POS") %>% 
  group_by(label, prediction) %>%
  tally() %>%
  mutate(pct = 100*(n / sum(n)))


