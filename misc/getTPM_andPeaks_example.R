# Get example data for testing different R-loop classification methods
library(tidyverse)

# Install RSeqR
remotes::install_github("Bishop-Laboratory/RSeqR", dependencies = "Imports", force = TRUE,
                         auth_token = "...")


# Annotate the RLoops
RLoops <- regioneR::toGRanges("misc/a423_rep1.broadPeak")
names(mcols(RLoops)) <- c("name", "score", "strand2", "pileup", "pVal", "qVal")

RLAnno <- RSeqR::annotateRLoops(RLoops)

# Get the TPM data
geneTPM <- readr::read_tsv("misc/a423_rep1.quant.sf") %>%
  select(Name, TPM)

save(geneTPM, RLAnno, file = "misc/geneTPM_and_RLAnno.rda", compress = "xz")






