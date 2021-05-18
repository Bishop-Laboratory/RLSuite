##### Correlation analysis #####
library(tidyverse)
library(ChIPpeakAnno)

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

#tagMatrix <- returnTagMatrix(int lower bound of conservation group, upper bound of conservation group, S4 peaks) 
returnTagMatrix <- function(lbound,ubound, rLoopRanges){
  #txdb<-hg38 required
  rangesLimited <- rLoopRanges[which(rLoopRanges$pct_cons >lbound & rLoopRanges$pct_cons <ubound )]
  
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  tagMatrixrangesLimited <- getTagMatrix(rangesLimited, windows = promoter)
  
  return(tagMatrixrangesLimited)
}

#Conservation groups are groups of R-loop windows within a specified range
#of percent score

#Creating tag maticies and metaplots for conservation groups in increments
#of 10% (not enough data to justify two groups for 80-90 and 90-100% conservation levels)

tagMatrix_rl_cons_0_10 <- returnTagMatrix(0,10,rl_cons)
tagMatrix_rl_cons_10_20 <- returnTagMatrix(10,20,rl_cons)
tagMatrix_rl_cons_20_30 <- returnTagMatrix(20,30,rl_cons)
tagMatrix_rl_cons_30_40 <- returnTagMatrix(30,40,rl_cons)
tagMatrix_rl_cons_40_50 <- returnTagMatrix(40,50,rl_cons)
tagMatrix_rl_cons_50_60 <- returnTagMatrix(50,60,rl_cons)
tagMatrix_rl_cons_60_70 <- returnTagMatrix(60,70,rl_cons)
tagMatrix_rl_cons_70_80 <- returnTagMatrix(70,80,rl_cons)
tagMatrix_rl_cons_80_100 <- returnTagMatrix(80,100,rl_cons)

plotAvgProf(tagMatrix_rl_cons_0_10, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_0_20")
plotAvgProf(tagMatrix_rl_cons_10_20, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_10_20")
plotAvgProf(tagMatrix_rl_cons_20_30, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_20_30")
plotAvgProf(tagMatrix_rl_cons_30_40, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_30_40")
plotAvgProf(tagMatrix_rl_cons_40_50, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_40_50")
plotAvgProf(tagMatrix_rl_cons_50_60, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_50_60")
plotAvgProf(tagMatrix_rl_cons_60_70, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_60_70")
plotAvgProf(tagMatrix_rl_cons_70_80, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_70_80")
plotAvgProf(tagMatrix_rl_cons_80_100, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3') rl_cons_80_100")

##Peak annitation and upsetplots

#The same 10% increments were used to group the windows

rl_cons_010_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >0 & rl_cons$pct_cons<10)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_1020_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >10 & rl_cons$pct_cons<20)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_2030_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >20 & rl_cons$pct_cons<30)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_3040_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >30 & rl_cons$pct_cons<40)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_4050_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >40 & rl_cons$pct_cons<50)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_5060_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >50 & rl_cons$pct_cons<60)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_6070_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >60 & rl_cons$pct_cons<70)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_7080_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >70 & rl_cons$pct_cons<80)], tssRegion=c(-3000, 3000),TxDb=txdb)
rl_cons_80100_anno <- annotatePeak(rl_cons[which(rl_cons$pct_cons >80 & rl_cons$pct_cons<100)], tssRegion=c(-3000, 3000),TxDb=txdb)

upsetplot(rl_cons_010_anno, vennpie = TRUE)
upsetplot(rl_cons_1020_anno, vennpie = TRUE)
upsetplot(rl_cons_2030_anno, vennpie = TRUE)
upsetplot(rl_cons_3040_anno, vennpie = TRUE)
upsetplot(rl_cons_4050_anno, vennpie = TRUE)
upsetplot(rl_cons_5060_anno, vennpie = TRUE)
upsetplot(rl_cons_6070_anno, vennpie = TRUE)
upsetplot(rl_cons_7080_anno, vennpie = TRUE)
upsetplot(rl_cons_80100_anno, vennpie = TRUE)

##Pathway Enrichment 

#Dotplots were commented out due to rendering issues
rl_cons_010_path <- enrichPathway(as.data.frame(rl_cons_010_anno)$geneId)
rl_cons_1020_path <- enrichPathway(as.data.frame(rl_cons_1020_anno)$geneId)
rl_cons_2030_path <- enrichPathway(as.data.frame(rl_cons_2030_anno)$geneId)
rl_cons_3040_path <- enrichPathway(as.data.frame(rl_cons_3040_anno)$geneId)
rl_cons_4050_path <- enrichPathway(as.data.frame(rl_cons_4050_anno)$geneId)
rl_cons_5060_path <- enrichPathway(as.data.frame(rl_cons_5060_anno)$geneId)
rl_cons_6070_path <- enrichPathway(as.data.frame(rl_cons_6070_anno)$geneId)
rl_cons_7080_path <- enrichPathway(as.data.frame(rl_cons_7080_anno)$geneId)
rl_cons_80100_path <- enrichPathway(as.data.frame(rl_cons_80100_anno)$geneId)

#dotplot(rl_cons_010_path,showCategory =20,title = "rl_cons 0-10% pctscore Pathways")
#dotplot(rl_cons_1020_path,showCategory =20,title = "rl_cons 10-20% pctscore Pathways")
#dotplot(rl_cons_2030_path,showCategory =20,title = "rl_cons 20-30% pctscore Pathways")
#dotplot(rl_cons_3040_path,showCategory =20,title = "rl_cons 30-40% pctscore Pathways")
#dotplot(rl_cons_4050_path,showCategory =20,title = "rl_cons 40-50% pctscore Pathways")
#dotplot(rl_cons_5060_path,showCategory =20,title = "rl_cons 50-60% pctscore Pathways")
#dotplot(rl_cons_6070_path,showCategory =20,title = "rl_cons 60-70% pctscore Pathways")
#dotplot(rl_cons_7080_path,showCategory =20,title = "rl_cons 70-80% pctscore Pathways")
#dotplot(rl_cons_80100_path,showCategory =20,title = "rl_cons 80-100% pctscore Pathways")

datatable(subset(as.data.frame(rl_cons_010_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 0-10 percent score"),
            filter = 'top')

datatable(subset(as.data.frame(rl_cons_1020_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 10-20 percent score"),
            filter = 'top')
datatable(subset(as.data.frame(rl_cons_2030_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 20-30 percent score"),
            filter = 'top')
datatable(subset(as.data.frame(rl_cons_3040_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 30-40 percent score"),
            filter = 'top')
datatable(subset(as.data.frame(rl_cons_4050_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 40-50 percent score"),
            filter = 'top')
datatable(subset(as.data.frame(rl_cons_5060_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 50-60 percent score"),
            filter = 'top')
datatable(subset(as.data.frame(rl_cons_6070_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 60-70 percent score"),
            filter = 'top')
datatable(subset(as.data.frame(rl_cons_7080_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 70-80 percent score"),
            filter = 'top')
datatable(subset(as.data.frame(rl_cons_80100_path),select = -c(geneID)),
            rownames = FALSE,
            caption = paste("pathway enrichment for R-loops: 80-100 percent score"),
            filter = 'top')

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

GT010 <-gene_table(0,10,rl_anno_2)
GT1020 <- gene_table(10,20,rl_anno_2)
GT2030 <- gene_table(20,30,rl_anno_2)
GT3040 <- gene_table(30,40,rl_anno_2)
GT4050 <- gene_table(40,50,rl_anno_2)
GT5060 <- gene_table(50,60,rl_anno_2)
GT6070 <- gene_table(60,70,rl_anno_2)
GT7080 <- gene_table(70,80,rl_anno_2)
GT80100 <- gene_table(80,100,rl_anno_2)

##Enrichment Analysis 

#Setting up enrichR
library(enrichR)
setEnrichrSite("Enrichr")
websiteLive<- TRUE

if(is.null(listEnrichrDbs()))
  print("DATABASES ARE NULL")
  websiteLive<-FALSE

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


