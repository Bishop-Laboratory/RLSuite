library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)


load("geneTPM_and_RLAnno.rda")
#imports geneTPM and RLAnno


# summing tpm to gene level  ----------------------------------------------

head(geneTPM)

geneTPM$TXID <- gsub("\\..*","", geneTPM$Name) 

sumtpm <- geneTPM %>% dplyr::select(c(TXID, TPM)) %>% group_by(TXID) %>% summarise(TPM=n())

head(sumtpm)

#don't need to summarise at the transcript level - should do this after mapping transcripts to genes 

# mapping transcripts to genes --------------------------------------------



anno <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = sumtpm$TXID, columns = c("TXID", "SYMBOL"),
                              keytype = "TXID")



df <- left_join(df, sumtpm) %>% dplyr::select(c(SYMBOL, TPM))

head(df)


# R loop wrangling --------------------------------------------------------

head(RLAnno)

rloops <- data.frame(SYMBOL = RLAnno$SYMBOL, RLOOP_PILEUP = RLAnno$pileup)

rloops <- rloops %>% group_by(SYMBOL) %>% summarise(RLoop_Pileup=n()) %>% drop_na()

df <- left_join(df, rloops)

head(df)

