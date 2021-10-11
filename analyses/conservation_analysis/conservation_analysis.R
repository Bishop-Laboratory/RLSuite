##### Correlation analysis #####

# Get the metadata
sample_metadata <- read_csv("misc/rmap_full_11_25.csv")

# Steps to get the data/hg38_10kb-widow_2.5kb-step_tiles.bed file ##

# Get the human genome
download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
              destfile = "misc/conservation_analysis/data/hg38.chrom.sizes")

# Divide the genome into 10kb windows with a 2.5kb step size
if (! file.exists("misc/conservation_analysis/data/hg38.10kb.tiles.bed")) {
  cmd <- paste0('bedtools makewindows -g misc/conservation_analysis/data/hg38.chrom.sizes ',
                '-w 10000 -s 2500 > misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.bed')
  system(cmd)
}

# Count number of overlaps with human experimental samples (removes control & non-human samples)
samples_to_intersect <- sample_metadata$sample_name[
  which(sample_metadata$genome == "hg38" &
          ! sample_metadata$Condition %in% c("WKKD", "IgG", "RNH", "ACTD"))
  ]

# Compile file paths for these peaks
peaks_to_intersect <- paste0('misc/conservation_analysis/data/final-peaks-unstranded/',
                             samples_to_intersect, '_hg38.unstranded.bed', collapse = " ")

# Count the overlaps between peaks and genome widows
cmd <- paste0('bedtools intersect -a misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.bed -b ',
              peaks_to_intersect,
              " -C > misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed")
if (! file.exists('misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed')) {
  system(cmd)
}

# Pivot the resulting table into a matrix
rl_cons_raw <- read_tsv("misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed",
                        col_names = c("seqnames", "start", "end", "intersected_file", "number_of_intersects"))
number_of_peak_files <- 252
rl_cons_mat <- rl_cons_raw %>%
  mutate(id = rep(paste0("window_", seq(length(rl_cons_raw$seqnames) / number_of_peak_files)), each = number_of_peak_files)) %>%
  mutate(number_of_intersects = ifelse(number_of_intersects > 0, 1, 0)) %>%
  mutate(intersected_file = paste0("peakfile_", intersected_file)) %>%
  select(id, intersected_file, number_of_intersects) %>%
  pivot_wider(names_from = intersected_file, values_from = number_of_intersects)
rl_cons_mat2 <- column_to_rownames(rl_cons_mat, var = "id")
rs <- rowSums(rl_cons_mat2)

window_key <- rl_cons_raw %>%
  mutate(id = rep(paste0("window_", seq(length(rl_cons_raw$seqnames) / number_of_peak_files)), each = number_of_peak_files)) %>%
  select(id, seqnames, start, end) %>%
  distinct(id, .keep_all = TRUE)
conservation_simple <- window_key %>%
  mutate(score = rs)
write_tsv(conservation_simple, path = "misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.conservation_levels.bed")

## Steps to do the conservation analysis -- start here ##

# Read in the results as a GRanges object
rl_cons <- read_tsv("misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.conservation_levels.bed", skip = 1,
                    col_names = c("name", "seqnames", "start", "end", "score"))
rl_cons <- toGRanges(as.data.frame(rl_cons))

# The score column of rl_cons represents the number of samples overlapping with that 10kb range
# See the histogram of scores
hist(rl_cons$score)
# Total number of samples is 252
total_number_of_samples <- 252

# Get the pct conservation
rl_cons$pct_cons <- 100 * (rl_cons$score / 252)

# See histogram of pct conservation
hist(rl_cons$pct_cons)

# Some windows have many R-loops overlapping with them, some have very few
top_ranges <- as.data.frame(rl_cons) %>%
  top_n(10, score) %>%
  rownames_to_column() %>%
  pull(rowname)
rl_cons[top_ranges,]  # Top 10 Ranges with the maximum score

# Continue with the rest of the analysis here...

##Dependencies 
#library(tidyverse)
#library(ChIPpeakAnno)
#library(ChIPseeker)
#library(ReactomePA)
#library(DT)

##Metaplot Analysis

#preparing annotations 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#metaplot <- metaplotViaTagMatrix(int vector lower bound of conservation group, int vector upper bound of conservation group, S4 peaks) 
metaplotViaTagMatrix <- function(lbound,ubound, rLoopRanges){
  #txdb<-hg38 required
  
  if (length(lbound)!= length(ubound)){
    print("Upper and lower bound vectors must be the same length!")
    return(NULL)
  }
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  
  for(i in 1:length(lbound)){
    rangesLimited <- rLoopRanges[which(rLoopRanges$pct_cons>lbound[i] & rLoopRanges$pct_cons <ubound[i])]
    tagMatrixrangesLimited <- getTagMatrix(rangesLimited, windows = promoter)
    label = paste(lbound[i],"-",ubound[i],"%",sep="")
    print(plotAvgProf(tagMatrixrangesLimited, xlim=c(-3000, 3000), xlab=label))
    
  }
  
 return(TRUE) 
}

lowerbound <- c(0,10,20,30,40,50,60,70,80)
upperbound <- c(10,20,30,40,50,60,70,80,100)
metaplotViaTagMatrix(lowerbound,upperbound,rl_cons)

##Peak annotation and upsetplots

#The same 10% increments were used to group the windows

#upsetplot <- annotateAndUpset(int vector lower bound of conservation group, int vector upper bound of conservation group, S4 peaks)
annotateAndUpset <- function(lbound,ubound,rLoopRanges){

    if (length(lbound)!= length(ubound)){
    print("Upper and lower bound vectors must be the same length!")
    return(NULL)
    }
    
    for(i in 1:length(lbound)){
      rl_annotated<- annotatePeak(rLoopRanges[which(rLoopRanges$pct_cons >lbound[i] &          rLoopRanges$pct_cons<ubound[i])], 
                                       tssRegion=c(-3000, 3000),TxDb=txdb)
      show(upsetplot(rl_annotated, vennpie = TRUE))
    }
    
}

annotateAndUpset(lowerbound,upperbound,rl_cons)

##Gene data tables 

#Each table contains the start, end, width, ENTREZID and description of the 
#genes found within each window.

annoData<-toGRanges(txdb)

rl_anno <- annotatePeakInBatch(rl_cons, 
                               AnnotationData=annoData, 
                               output="overlapping",
                               bindingRegion=c(-3000, 3000))

rl_anno_2 <- addGeneIDs(rl_anno,
                        "org.Hs.eg.db",
                        IDs2Add = c("symbol","genename"),
                        feature_id_type = "entrez_id")

names(rl_anno_2) = make.names(names(rl_anno_2),unique = TRUE)

#datatable <- gene_table(int lower bound of conservation group, int upper bound of conservation group, S4 annotated peaks with gene Ids)
gene_table<- function(lbound,ubound, rLoopRanges){
  rl_limited <- rLoopRanges[which(rLoopRanges$pct_cons >lbound & rLoopRanges$pct_cons <ubound )]
  return(
    datatable(subset(as.data.frame(rl_limited),select = -c(width,strand,score,peak)),
          caption = paste(lbound,"-",ubound,"% score",sep= ""),
          extensions = 'FixedColumns',
          options = list(
          dom = 't',
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 2, rightColumns = 2),
          scrollCollapse = TRUE))
  )
}

if(length(upperbound) == length(lowerbound)){
  for(i in 1:length(upperbound)){
    show(gene_table(lowerbound[i], upperbound[i], rl_anno_2))
  }
}

##Enrichment Analysis 

#Setting up enrichR

#ChEA database was used to create datables 

#80-100%
rl80100 <- rl_anno_2[which(rl_anno_2$pct_cons >80 & rl_anno_2$pct_cons <100 )]
rl80100_symbols <- unique(rl80100$symbol)
enriched80100 <- enrichr(rl80100_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e80100<-enriched80100[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e80100[e80100$P.value>=0.01 & e80100$P.value<=0.05,],
          caption = "80-100% conservation",
          options = list(
          dom = 't',
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1),
          scrollCollapse = TRUE))

#70-80%
rl7080 <- rl_anno_2[which(rl_anno_2$pct_cons >70 & rl_anno_2$pct_cons <80 )]
rl7080_symbols <- unique(rl7080$symbol)
enriched7080 <- enrichr(rl7080_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e7080<-enriched7080[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e7080[e7080$P.value>=0.01 & e7080$P.value<=0.05,],
          options = list(
          dom = 't',
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1),
          scrollCollapse = TRUE))

#60-70%
rl6070 <- rl_anno_2[which(rl_anno_2$pct_cons >60 & rl_anno_2$pct_cons <70 )]
rl6070_symbols <- unique(rl6070$symbol)
enriched6070 <- enrichr(rl6070_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e6070<-enriched6070[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e6070[e6070$P.value>=0.01 & e6070$P.value<=0.05,],
          options = list(
          dom = 't',
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1),
          scrollCollapse = TRUE))

#50-60%
rl5060 <- rl_anno_2[which(rl_anno_2$pct_cons >50 & rl_anno_2$pct_cons <60 )]
rl5060_symbols <- unique(rl5060$symbol)
enriched5060 <- enrichr(rl5060_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e5060<-enriched5060[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e5060[e5060$P.value>=0.01 & e5060$P.value<=0.05,],
          options = list(
          dom = 't',
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1),
          scrollCollapse = TRUE))

#40-50%
rl4050 <- rl_anno_2[which(rl_anno_2$pct_cons >40 & rl_anno_2$pct_cons <50 )]
rl4050_symbols <- unique(rl4050$symbol)
enriched4050 <- enrichr(rl4050_symbols, c("ChEA_2016"))
e4050<-enriched4050[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e4050[e7080$P.value>=0.01 & e7080$P.value<=0.05,],
          options = list(
          dom = 't',
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1),
          scrollCollapse = TRUE))

#30-40%
rl3040 <- rl_anno_2[which(rl_anno_2$pct_cons >30 & rl_anno_2$pct_cons <40 )]
rl3040_symbols <- unique(rl3040$symbol)
enriched3040 <- enrichr(rl3040_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e3040<-enriched3040[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e3040[e3040$P.value>=0.01 & e3040$P.value<=0.05,],
          options = list(
          dom = 't',
          scrollX = TRUE,
          fixedColumns = list(leftColumns = 1),
          scrollCollapse = TRUE))

#20-30%
rl2030 <- rl_anno_2[which(rl_anno_2$pct_cons >20 & rl_anno_2$pct_cons <30 )]
rl2030_symbols <- unique(rl2030$symbol)
enriched2030 <- enrichr(rl2030_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e2030<-enriched2030[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e2030[e2030$P.value>=0.01 & e2030$P.value<=0.05,])

#10-20%
rl1020 <- rl_anno_2[which(rl_anno_2$pct_cons >10 & rl_anno_2$pct_cons <20 )]
rl1020_symbols <- unique(rl1020$symbol)
enriched1020 <- enrichr(rl1020_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e1020<-enriched1020[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e1020[e1020$P.value>=0.01 & e1020$P.value<=0.05,])

#<=10%
rl10 <- rl_anno_2[which(rl_anno_2$pct_cons <10 )]
rl10_symbols <- unique(rl10$symbol)
enriched10 <- enrichr(rl10_symbols, c("ChEA_2016","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018"))
e10<-enriched10[["ChEA_2016"]]
#datatable(e[e$P.value<=0.05,])
datatable(e10[e10$P.value>=0.01 & e10$P.value<=0.05,])


