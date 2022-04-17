library(tidyverse)
library(enrichR)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
resloc <- "results/New_Figres/"
num_sel <- 6
chain <- "data/hg38_to_hg19.chain.gz"
chn <- rtracklayer::import.chain(gsub(chain, pattern = "\\.gz", replacement = ""))

## hg19 windows
gen <- valr::read_genome("tmp/hg19.chrom.sizes")
gen %>% mutate(start = 0, end = size) -> dd

## Download hg19 tracks needed for visualizations
bwlist <- list(
  "iPSCs" = list(
    "ATAC" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136857/suppl/GSE136857_iPSC_IDR.bw",
    "EZH2" = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1180nnn/GSM1180135/suppl/GSM1180135_AM-iPS-6.EZH2.v2m4bhg19ren.bg.bw",
    "CTCF" = "https://www.encodeproject.org/files/ENCFF832IJN/@@download/ENCFF832IJN.bigWig",
    "GRO-Seq" = "https://rlbase-data.s3.amazonaws.com/misc/iPSCs_GROSeq.bw",
    "H3K27ac" = "https://www.encodeproject.org/files/ENCFF833FFJ/@@download/ENCFF833FFJ.bigWig",
    "POLR2A" = "https://www.encodeproject.org/files/ENCFF919VDD/@@download/ENCFF919VDD.bigWig",
    "H3K27me3" = "https://www.encodeproject.org/files/ENCFF645JYS/@@download/ENCFF645JYS.bigWig",
    "H3K4me1" = "https://www.encodeproject.org/files/ENCFF415ZUP/@@download/ENCFF415ZUP.bigWig",
    "MapR" = "https://rlbase-data.s3.amazonaws.com/misc/iPSC_MapR_comb.bw",
    "S96-RLRegions" = "https://rlbase-data.s3.amazonaws.com/misc/rlregions_S96.hg19.bw.bw",
    "dRNH-RLRegions" = "https://rlbase-data.s3.amazonaws.com/misc/rlregions_dRNH.hg19.bw.bw"
  )
)
lapply(names(bwlist), function(cell) {
  message(cell)
  dir.create(file.path("tmp", cell), showWarnings = FALSE)
  lapply(
    seq(bwlist[[cell]]), function(i) {
      message(i)
      dest <- file.path(
        "tmp", cell, paste0(cell, "__", names(bwlist[[cell]])[i], ".bw")
      )
      if (! file.exists(dest)) {
        download.file(
          bwlist[[cell]][[i]], 
          destfile = dest
        )
      }
    }
  )
})

## Get enhancers from ChromHMM
fls <- list(
  "HEK293" = "https://www.encodeproject.org/files/ENCFF476TTU/@@download/ENCFF476TTU.bed.gz",
  "HCT116" = "https://www.encodeproject.org/files/ENCFF513PJK/@@download/ENCFF513PJK.bed.gz",
  "iPSCs" = "https://www.encodeproject.org/files/ENCFF115RIR/@@download/ENCFF115RIR.bed.gz"
)
flshg19 <- list(
  "HEK293" = "https://www.encodeproject.org/files/ENCFF402EKL/@@download/ENCFF402EKL.bed.gz",
  "iPSCs" = "https://www.encodeproject.org/files/ENCFF197OWG/@@download/ENCFF197OWG.bed.gz"
)

pkschrom <- lapply(
  fls, read_tsv,
  col_names = c(
    "chrom", "start", "end", "name",
    "score", "strand", "start2", "end2", "color"
  )
)
pkschrom <- lapply(
  pkschrom, function(x) {
    # Extract chromhmm predicted enhancers
    filter(x, str_detect(name, "Enh")) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  }
)


## Get GH enhancers
gh <- readRDS("data/genehancer_processed.rds")
ghfull <- gh %>%
  unnest(cols = data)
## Find the distal enhancers
ghgr <- gh %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
ghanno <- ghgr %>%
  ChIPseeker::annotatePeak(TxDb = txdb, tssRegion = c(-3000, 3000), verbose = T)
dEnh <- ghgr[which(ghanno@detailGenomicAnnotation$Intergenic), ]

###$$$### ANalyze iPSCs

ct <- "iPSCs"

### GRO-Seq
infile <- paste0("analyses/extension_of_dRNH_S96_analysis/", ct, "_GROSeq.bdg")
infile2 <- paste0("analyses/extension_of_dRNH_S96_analysis/", ct, "_GROSeq.sort.bdg")
outfile <- paste0("analyses/extension_of_dRNH_S96_analysis/", ct, "_GROSeq.bw")
if (! file.exists(outfile)) {
  
  # Make genome windows
  wins <- valr::bed_makewindows(x = dd, win_size = 50)
  
  # Read in GRO-Seq
  gro1 <- valr::read_bedgraph(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3271nnn/GSM3271001/suppl/GSM3271001_iPSC_11104_A.bedGraph.gz"
  )
  gro2 <- valr::read_bedgraph(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3271nnn/GSM3271002/suppl/GSM3271002_iPSC_11104_B.bedGraph.gz"
  )
  
  # Map to windows
  gros <- parallel::mclapply(
    list(gro1, gro2), function(x) {
      valr::bed_map(wins, x, score2 = mean(value)) %>% 
        mutate(score2 = ifelse(is.na(score2), 0, score2))
    }, mc.cores = 2
  )
  
  # Average scores
  GRO <- select(gros[[1]], chrom, start, end)
  GRO$score <- (gros[[1]]$score2 + gros[[2]]$score2) / 2
  
  # Write to disk and convert to bigWig
  GRO %>% write_tsv(infile, col_names = FALSE)
  system(paste0("tail -n +2 ", infile, " | sort -k1,1 -k2,2n > ", infile2))
  system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedGraphToBigWig ", infile2, " tmp/hg19.chrom.sizes ", outfile))
  aws.s3::put_object(file = outfile, object = paste0("misc/", ct, "_GROSeq.bw"), bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
  
}

### Analyze BiV enhancers

# Get BiV Enhancer annotations from chromHMM
chromenh <- pkschrom[[ct]]
chrombiv <- chromenh[chromenh$name == "EnhBiv",]

## Load in the MapR data
dat <- RLHub::rlbase_samples()
ctl <- dat %>%
  filter(
    ip_type == "dRNH",
    genome == "hg38",
    mode == "MapR",
    tissue == ct,
    label == "POS",
    prediction == "POS",
    numPeaks > 5000
  ) %>%
  slice_sample(n=5)
grs <- lapply(ctl$rlsample, function(x) {
  GenomicRanges::GRanges(RLSeq::RLRangesFromRLBase(x))
})
names(grs) <- ctl$rlsample


## Analyze feature distribution
pal <- lapply(
  grs, ChIPseeker::annotatePeak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000), verbose = T
)

## Get intergenic MapR peaks
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)
gri <- GenomicRanges::reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]


## Distance to bivalent enhancer
pltd <- ChIPpeakAnno::binOverFeature(
  gris[[1]], gris[[2]], gris[[3]], gris[[4]],
  annotationData = chrombiv,
  radius = 5000, nbins = 100, FUN = length, errFun = 0,
  xlab = "Distance from Bivalent Enhancer (bp)", ylab = "count"
)
plt <- pltd %>%
  as.data.frame() %>%
  rownames_to_column("pos") %>%
  as_tibble() %>%
  mutate(pos = as.numeric(pos)) %>%
  dplyr::rename(
    "SRX10505690" = 2, 
    "SRX10505691" = 3, 
    "SRX10505696" = 4,
    "SRX10505697" = 5
  ) %>%
  pivot_longer(cols = !contains("pos")) %>%
  ggplot(aes(x = pos, y = value, color = name)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .25) +
  geom_point() +
  geom_line() +
  ggtitle("Peak pileup around bivalent enhancers", subtitle = "iPSCs MapR (dRNH) peaks") +
  ylab("Peak density") +
  xlab("Distance to bivalent enhancer (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2")
plt
ggsave(plt, filename = file.path(resloc, paste0(ct, "_pileup_bivEnh.png")))


## Find the specific Enh
# Get the intergenic peaks from both iPSC samples & find overlap
pki <- pkschrom[[ct]] %>%
  ChIPseeker::annotatePeak(TxDb = txdb, tssRegion = c(-3000, 3000), verbose = TRUE)
pki <- pkschrom[[ct]][pki@detailGenomicAnnotation$Intergenic, ]
pki <- pki[!pki$name %in% c("EnhG1", "EnhG2"), ]

# Transfer labels to the dEnh from GH
olpe <- ChIPpeakAnno::findOverlapsOfPeaks(
  pki, dEnh[dEnh$score > .5, ]
)
colnames(olpe$overlappingPeaks$`pki///dEnh`)[7] <- "pName"
dEnh2 <- olpe$overlappingPeaks$`pki///dEnh` %>%
  select(name, pName) %>%
  as_tibble() %>%
  group_by(name) %>%
  mutate(n = n()) %>%
  filter(n == 1) %>%
  select(-n) %>%
  inner_join(as.data.frame(dEnh)) %>%
  as.data.frame() %>%
  ChIPpeakAnno::toGRanges()

# Analyze dEnh and Peaks
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh2, gri)
d1 <- olde3$overlappingPeaks$`dEnh2///gri`$pName %>%
  tibble() %>%
  dplyr::rename(enh = 1) %>%
  group_by(enh) %>%
  tally() %>%
  mutate(`dRNH-accessible` = n / sum(n)) %>%
  select(-n)
plt <- olde3$all.peaks$dEnh2$pName %>%
  tibble() %>%
  dplyr::rename(enh = 1) %>%
  group_by(enh) %>%
  tally() %>%
  mutate(`Total dEnh` = n / sum(n)) %>%
  select(-n) %>%
  inner_join(d1) %>%
  pivot_longer(cols = !matches("^enh")) %>%
  ggplot(aes(fill = enh, y = value, x = name)) +
  geom_col(color="black") +
  coord_flip() +
  xlab(NULL) +
  ylab("Proportion of distal enhancers") +
  ggtitle(
    paste0(ct, " distal enhancer profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  scale_fill_manual(
    values = c(
      "EnhA1" = "#dbdbdb",
      "EnhA2" = "#ededed",
      "EnhBiv" = "#4141bf",
      "EnhWk" = "#c9c7c7"
    ), labels = c(
      "Active (1)", 
      "Active (2)", 
      "Bivalent",
      "Weak"
    )
  ) +
  theme_bw(base_size = 14) +
  guides(fill = guide_legend(title = " ")) 
plt
ggsave(plt, filename = file.path(resloc, paste0(ct, "_dENH_types_prop.png")))

# For these overlapping peaks, what are they?
num_sel <- 6
dENH_gri <- olde3$overlappingPeaks$`dEnh2///gri`
ghfull %>%
  filter(
    name %in% dENH_gri$peaks1[dENH_gri$pName %in% c("EnhBiv")],
    gene_scores > 10
  ) -> dd
unique(dd$genes) -> genesNow
genesNow
eres <- enrichr(genesNow, databases = "ChEA_2016")
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>%
  filter(row_number() <= num_sel) %>%
  pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms) %>%
  select(Term, combined_score = Combined.Score, padj = Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>%
  mutate(
    Term = factor(Term,
                  levels = unique(Term)
    )
  ) %>%
  filter(!is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col(fill="#F58776") +
  coord_flip() +
  theme_bw(base_size = 14) +
  xlab(NULL) +
  ylab("Enrichment (Combined Score)") +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in ChEA Database")
plt
ggsave(plt, filename = file.path(resloc, paste0(ct, "_dRNH-dEnh_Biv_chea.png")))

# CellMarker
num_sel <- 6
dENH_gri <- olde3$overlappingPeaks$`dEnh2///gri`
ghfull %>%
  filter(
    name %in% dENH_gri$peaks1[dENH_gri$pName %in% c("EnhBiv")],
    gene_scores > 10
  ) -> dd
unique(dd$genes) -> genesNow
genesNow
eres <- enrichr(genesNow, databases = "CellMarker_Augmented_2021")
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>%
  filter(row_number() <= num_sel) %>%
  pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms) %>%
  select(Term, combined_score = Combined.Score, padj = Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>%
  mutate(
    Term = factor(Term,
                  levels = unique(Term)
    )
  ) %>%
  filter(!is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col(fill="#F58776") +
  coord_flip() +
  theme_bw(base_size = 14) +
  xlab(NULL) +
  ylab("Enrichment (Combined Score)") +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in ChEA Database")
plt
ggsave(plt, filename = file.path(resloc, paste0(ct, "_dRNH-dEnh_Biv_chea.png")))

eres <- enrichr(genesNow, databases = "ARCHS4_TFs_Coexp")
num_sel <- 12
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>%
  filter(row_number() <= num_sel) %>%
  pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms, P.value < 0.01) %>%
  select(Term, combined_score = Combined.Score, padj = Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>%
  mutate(
    Term = gsub(Term, pattern = " human tf ARCHS4 coexpression", replacement = ""),
    Term = factor(Term,
                  levels = unique(Term)
    )
  ) %>%
  filter(!is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col(fill="#F58776") +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab("Enrichment (Combined Score)") +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in ARCHS4 TF Database")
plt
ggsave(plt, filename = file.path(resloc, paste0(ct, "_dRNH-dEnh_Biv_archs4.png")))


#### Analyze bivalent enhancer clusters

lenmin <- 1

# Calculate clusters for all non-bivalent enhancers
chromenh <- pkschrom[[ct]] 
chromenh <- chromenh[chromenh$name != "EnhBiv", ]
chromenhhg19 <- rtracklayer::liftOver(chromenh, chain = chn) %>% unlist()
chromenhhg19 <- unique(chromenhhg19)
bigchromenh <- (chromenhhg19 + 5000) %>% GenomicRanges::reduce(with.revmap=T)
len <- sapply(bigchromenh$revmap, length)
bigchromenh[len > lenmin, ] %>%  #[sample(length(bigchromenh), 10000),] %>%
  rtracklayer::export(
    con = paste0("tmp/enh_clusters__", ct, ".bed")
  )

# Calculate clusters for bivalent enhancers
chromenh <- pkschrom[[ct]] 
chromenhbiv <- chromenh[chromenh$name == "EnhBiv", ]
chromenhbivhg19 <- rtracklayer::liftOver(chromenhbiv, chain = chn) %>% unlist()
chromenhbivhg19 <- unique(chromenhbivhg19)
bigchromenhbiv <- (chromenhbivhg19 + 5000) %>% GenomicRanges::reduce(with.revmap=T)
len <- sapply(bigchromenhbiv$revmap, length)
bigchromenhbiv[len > lenmin, ] %>%  
  rtracklayer::export(
    con = paste0("tmp/enh_biv_clusters__", ct, ".bed")
  )

# Get track for hg19 chromHMM
infile <- paste0("tmp/", ct, "_chromhmm_enhbiv.bed")
infile2 <- paste0("tmp/", ct, "_chromhmm_enhbiv.sort.bed")
outfile <- paste0("tmp/", ct, "_chromhmm_enhbiv.sort.bb")
if (! file.exists(outfile)) {
  hg19_chromhmm <- valr::read_bed12(flshg19[[ct]])
  hg19_chromhmm %>% filter(name == "EnhBiv") %>%
    write_tsv(infile, col_names = F)
  system(paste0("tail -n +2 ", infile, " | sort -k1,1 -k2,2n > ", infile2))
  system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedToBigBed ", infile2, " tmp/hg19.chrom.sizes ", outfile))
  aws.s3::put_object(file = outfile, object = paste0("misc/", ct, "_hg19_enhbiv.bb"), bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
}

## Run computeMatrix
bwsct <- lapply(seq(bwlist[[ct]]), function(i) {
  file.path(
    "tmp", ct, paste0(ct, "__", names(bwlist[[ct]])[i], ".bw")
  )
})
names(bwsct) <- names(bwlist[[ct]])
outmat <- paste0('tmp/enh_clusters__', ct, '.mat.gz')
sampleOrder <- c(
  "GRO-Seq", 
  "POLR2A",
  "MapR", 
  "dRNH-RLRegions", 
  "S96-RLRegions",
  "CTCF",
  "EZH2",
  "ATAC",
  "H3K27me3",
  "H3K4me1",
  "H3K27ac"
)
system(
  paste0(
    '~/miniconda3/condabin/conda run -n rlbaseData computeMatrix scale-regions ',
    '-S ',
    paste(
      sapply(sampleOrder, function(x) {
        bwsct[[x]]
      }), collapse = " "
    ),
    ' -R ', paste0("tmp/enh_biv_clusters__", ct, ".bed"), " ", paste0("tmp/enh_clusters__", ct, ".bed"),
    ' -m 15000 -a 15000 -b 15000 -bs 500 -p 44 --skipZeros -o ',
    outmat
  )
)

## Parse matrix output
n <- length(sampleOrder)

# Read in matrix & metadata
mat <- read_tsv(outmat, skip = 1, col_names = FALSE)
meta <- read_lines(outmat)[1]
meta <- gsub(meta, pattern = "@", replacement = "")
meta <- jsonlite::parse_json(meta)

# Get the column names from IP name and bp position
is <- seq(meta$sample_boundaries)
cnames <- lapply(
  is[-length(is)], function(i) {
    names(bwsct)[i]
    poss <- (1+meta$sample_boundaries[[i]]):(meta$sample_boundaries[[i+1]])
    paste0(names(bwsct)[i], "__", poss)
  }
) %>% unlist()
colnames(mat)[1:6] <- c("chrom", "start", "end", "name", "score", "strand")
colnames(mat)[7:length(colnames(mat))] <- cnames

# Convert to long-form and extract names, position & handle NaNs
mlong <- mat %>% pivot_longer(cols = contains("__"), names_to = "NAMES")
lpat <- "(.+)__(.+)"
mlong <- mlong %>% 
  mutate(
    ip= str_replace(NAMES, pattern = lpat, replacement = "\\1"),
    pos = str_replace(NAMES, pattern = lpat, replacement = "\\2"),
    value = as.numeric(case_when(value == "NaN" ~ 0, TRUE ~ value))
  )

# Get the library size for normalization
sums <- mlong %>%
  group_by(ip) %>% 
  summarise(val=sum(value))

# Library-size norm
mm <- function(x) (x-min(x))/(max(x)-min(x))
mlong2 <- inner_join(mlong, sums)
mlong3 <- mlong2 %>%
  group_by(ip) %>% 
  mutate(
    value_scale = round(1000000*(value / val), digits = 6),
    # value_scale = round((100*as.numeric(scale(value_scale, center = T))), digits=6)
    # value_scale = as.numeric(scale(value, center = T)),
    # value_scale = round((1000*mm(value)), digits = 6)
  )

# Reformat to wide form data
mwide <- mlong3 %>% 
  ungroup() %>% 
  select(-value, -ip, -pos, -val) %>% 
  pivot_wider(names_from = NAMES, values_from = value_scale)

# Save and convert back to deepTools format
outmat2 <- paste0('tmp/enh_clusters__', ct, '.2.mat')
outmat3 <- paste0('tmp/enh_clusters__', ct, '.3.mat')
write_tsv(mwide, file = outmat2, col_names = FALSE)
system(
  paste0(
    "zcat ", outmat, " | head -n 1 > ", outmat3, " && cat ", outmat2, 
    " >> ", outmat3, " && gzip -f ", outmat3
  )
)

# Get tornado plot
outres <- file.path(resloc, paste0(ct, "__biv_tornado.pdf"))
system(
  (
    paste0(
      '~/miniconda3/condabin/conda run -n rlbaseData plotHeatmap --matrixFile ',
      paste0(outmat3, ".gz"), ' --outFileName ', outres,
      ' --samplesLabel "', paste0(sampleOrder, collapse = '" "'), '"',
      " -z Bivalent All ",
      " --zMax 2.2 --zMin 0.5 --yMax 2.2 --yMin 0.2 ",
      '--startLabel	"start" --endLabel "end" --heatmapWidth 4.2 --heatmapHeight 10 --colorMap "inferno"'
    )
  )
)
print("DONEE")



