library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)


load("geneTPM_and_RLAnno.rda")
#imports geneTPM and RLAnno data


# summarising then mapping transcripts to genes --------------------------------------------

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

#create rloops df from example data using SYMBOL and pileup 
rloops <- data.frame(SYMBOL = RLAnno$SYMBOL, RLOOP_PILEUP = RLAnno$pileup)

#groupby SYMBOL & summarise, returning mean of pileup values linked to each symbol) 
rloops <- rloops %>% group_by(SYMBOL) %>% summarise(RLoop_Pileup_Mean=mean(RLOOP_PILEUP)) 

#merge rloops with df, creating df with symbol, tpm & pileup mean columns 
df <- left_join(df, rloops)


