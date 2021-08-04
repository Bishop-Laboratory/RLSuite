setwd("C:\\Users\\annaj\\OneDrive - University of Bristol\\Extracurricular\\UTSA")

library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)


load("geneTPM_and_RLAnno.rda")
#imports geneTPM and RLAnno data


# summarising then mapping transcripts to genes --------------------------------------------

head(geneTPM)

#put names into standard TXID format
geneTPM$TXID <- gsub("\\..*","", geneTPM$Name) 

#create df of TXID & TPM
df <- geneTPM %>% dplyr::select(c(TXID, TPM))

#find corresponding symbols for TXIDs
anno <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = df$TXID, columns = c("TXID", "SYMBOL"),
                              keytype = "TXID")

#create dataframe of symbols and TPM, summarising to gene level
df <- as.data.frame(left_join(df, anno)) %>% dplyr::select(c(SYMBOL, TPM)) %>% group_by(SYMBOL) %>% summarise(TPM = sum(TPM))



# R loop wrangling --------------------------------------------------------

head(RLAnno)

rloops <- data.frame(SYMBOL = RLAnno$SYMBOL, RLOOP_PILEUP = RLAnno$pileup)

rloops <- rloops %>% group_by(SYMBOL) %>% summarise(RLoop_Pileup=n()) %>% drop_na()

df <- left_join(df, rloops)

head(df)

df2['TPM','ABCA2']
