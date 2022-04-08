## This is an extension of the analysis from the original paper 
## Purpose is to characterize S96 and dRNH peaks and dig deeper into them

library(tidyverse)

#### MAGIC ####

ndrnh <- 56  # Number of dRNH passing QC filters
ns96 <- 191  # Number of S9.6 passing QC filters
annFull <- RLHub::annots_full_hg38()

###############

##### Analysis of peak sizes and locations dRNH vs S9.6 #####

#### Step #1: Find the difference in sizes between shared/non-shared consensus peaks

## Read in peaks
drnh_cons <- read_tsv(
  file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_dRNH.narrowPeak"),
  col_names = c(
    "chrom", "start", "end", "name", "score", 
    "strand", "signalVal", "pVal", "qVal", "peak"
  ), 
  show_col_types = FALSE, progress = FALSE
) 
s96_cons <- read_tsv(
  file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_S96.narrowPeak"),
  col_names = c(
    "chrom", "start", "end", "name", "score",
    "strand", "signalVal", "pVal", "qVal", "peak"
  ), 
  show_col_types = FALSE, progress = FALSE
)
## Read in genome info
genome <- read_tsv(
  "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", 
  col_names = c("chrom", "size"), show_col_types = FALSE, progress = FALSE
)
genes <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, 
  keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86), 
  columns = c("SEQNAME", "GENESEQSTART", "GENESEQEND", "SYMBOL")
) %>% dplyr::rename(chrom = SEQNAME, start = GENESEQSTART, end = GENESEQEND) %>%
  filter(! grepl(GENEID, pattern = "LRG.+")) %>%
  select(-GENEID) %>% mutate(chrom = paste0("chr", chrom)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

## Get the dRNH-only and S96-only along with shared
# TODO: Re-check this using ChIP peak anno
# int <- valr::bed_intersect(s96_cons, drnh_cons)
# olol <- ChIPpeakAnno::findOverlapsOfPeaks(GenomicRanges::makeGRangesFromDataFrame(s96_cons), GenomicRanges::makeGRangesFromDataFrame(drnh_cons))
# ChIPpeakAnno::makeVennDiagram(olol)

cons2 <- list(
  "dRNH-only" = drnh_cons[! drnh_cons$name %in% int$name.y,],
  "S9.6-only" = s96_cons[! s96_cons$name %in% int$name.x,],
  "dRNH-shared" = drnh_cons[drnh_cons$name %in% int$name.y,],
  "S9.6-shared" = s96_cons[s96_cons$name %in% int$name.x,],
  "S9.6" = s96_cons
)


## Density plot of peak sizes
minmax <- function(x) (x- min(x)) /(max(x)-min(x))
pltdat <- bind_rows(cons2, .id = "group") %>%
  mutate(
    width = end - start
  ) %>% 
  mutate(
    pct_cons = case_when(
      grepl(group, pattern = "dRNH") ~ 100 * ((score / 10) / ndrnh),
      grepl(group, pattern = "S9.6") ~ 100 * ((score / 10) / ns96)
    ),
    pct_rank = minmax(rank(pct_cons)),
    group2 = ifelse(grepl(group, pattern = "S9.6"), "S9.6", group)
  )
mu <- pltdat %>% group_by(group2) %>% summarise(mu=mean(width))
mu
plt <- ggplot(pltdat, aes(y = width, x=group2, fill = group2)) +
  geom_violin(width = .6) +
  geom_boxplot(width = .2) +
  # ggridges::geom_density_ridges(alpha=.3) +
  scale_y_log10() +
  theme_classic(base_size = 18) +
  ylab("Peak width (bp)") +
  xlab(NULL) +
  ggtitle("Peak width across groups", subtitle = "Consensus peaksets") +
  guides(fill=guide_legend(title=NULL), color=guide_legend(title=NULL)) +
  scale_fill_manual(values = c("dRNH-only" = "#e2a3e3", "dRNH-shared" = "#ccb6cc", "S9.6" = "#82d0e8")) +
  ggpubr::stat_compare_means(
    comparisons = list(c("dRNH-only", "dRNH-shared")),
    method.args = list(alternative="less")
  ) +
  ggpubr::rremove("legend")
plt

## Density plot of peak strengths
mu <- pltdat %>% group_by(group2) %>% summarise(mu=mean(pct_cons))
mu
plt <- ggplot(pltdat, aes(y = pct_cons, x=group2, color = group2)) +
  geom_violin(width = .6) +
  geom_boxplot(width = .2) +
  # ggridges::geom_density_ridges(alpha=.3) +
  scale_y_log10() +
  theme_classic(base_size = 18) +
  ylab("Peak width (bp)") +
  xlab(NULL) +
  ggtitle("Peak conservation % across groups", subtitle = "Consensus peaksets") +
  guides(fill=guide_legend(title=NULL), color=guide_legend(title=NULL)) +
  ggpubr::stat_compare_means(
    comparisons = list(c("dRNH-only", "dRNH-shared")),
    method.args = list(alternative="less")
  )
plt


#### Step #2: Proportion test
## What prop of dRNH-detected R-loops are also detected by S9.6?
## What prop of S9.6-detected R-loops are also detected by dRNH?
## Check with pie charts or Bar charts
pltdat %>% 
  filter(group != "S9.6") %>% 
  group_by(group) %>% 
  tally() %>% 
  mutate(group2 = gsub(group, pattern = "\\-.+", replacement = ""),
         group3 = gsub(group, pattern = ".+\\-", replacement = "")) %>% 
  group_by(group2) %>% 
  mutate(prop = n / sum(n)) %>% 
  ggplot(aes(x = group2, y = prop, fill=group3)) +
  geom_col() +
  coord_flip() +
  theme_gray(base_size = 14) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  ggtitle("Proportion of S9.6 and dRNH peaks found by both IPs vs not")



#### Step #3: Analysis of Tx Features
### What are the features overlapping consensus peaks in each group?
### How does this relate to peak size/strength?

cons3 <- cons2[! names(cons2) %in% c("S9.6-only", "S9.6-shared")]

## Compare samples around genomic features
consInt2 <- lapply(cons3, function(x) {
  # Call summits
  x <- x %>%
    dplyr::rename(so = start, 
                  eo = end) %>%
    mutate(
      start = so + peak - 250,
      end = so + peak + 250
    )
  # Annotate genes
  xint <- x %>%
    relocate(start, .after = chrom) %>%
    relocate(end, .after = chrom) %>%
    select(-so, -eo) %>%
    valr::bed_intersect(y = genes) %>%
    select(name = name.x, SYMBOL = SYMBOL.y) %>%
    right_join(x, by = "name") 
  return(xint)
}) %>% bind_rows(.id = "group")

# Annotate TTS, TSS, and Gene Body
txfeat_cols <- setNames(
  rev(
    c(
      "#b36262", "#c4cc6a", "#d17bdb",
      "#4bc6d6", "#83d647", "#9494e3", "#7d7d7d"
    )
  ),
  nm = rev(
    c(
      "TSS", "fiveUTR", "Exon",
      "Intron", "threeUTR", "TTS", "Intergenic"
    )
  )
)
txfeats <- annFull[sapply(
  names(annFull), 
  function(x) grepl(x, pattern = "Transcript_Features.*")
)]
txfeats <- lapply(names(txfeats), function(txname) {
  x <- txfeats[[txname]]
  x$db <- gsub(txname, pattern = "(.+)__(.+)", replacement = "\\1")
  x$type <- gsub(txname, pattern = "(.+)__(.+)", replacement = "\\2")
  return(x)
}) %>% bind_rows()

# With shared
oltx2 <- consInt2 %>%
  mutate(source = group) %>%
  select(chrom, start, end, name, source) %>%
  valr::bed_intersect(txfeats)
oltxsum2 <- oltx2 %>%
  group_by(name.x) %>%
  summarise(
    type = ifelse(
      "TSS" %in% type.y, "TSS", 
      ifelse("TTS" %in% type.y, "TTS",
             ifelse("fiveUTR" %in% type.y, "fiveUTR", 
                    ifelse("threeUTR" %in% type.y, "threeUTR", 
                           ifelse("Exon" %in% type.y, "Exon", "Intron"))))
    )
  )

oltxres2 <- oltxsum2 %>%
  dplyr::rename(name = name.x) %>%
  full_join(consInt2, by = "name") %>%
  mutate(source = group) %>%
  select(name, source, type) %>%
  mutate(type = ifelse(is.na(type), "Intergenic", type)) %>%
  distinct(name, .keep_all = TRUE)
plt <- oltxres2 %>%
  group_by(source, type) %>%
  tally() %>%
  mutate(
    n = 100*n / sum(n),
    type = factor(type, levels = rev(c(
      "TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "Intergenic"
    )))
  )  %>%
  filter(source != "Shared") %>% 
  mutate(source = factor(source, levels = rev(unique(.$source)))) %>% 
  ggplot(aes(x = source, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  coord_flip() +
  ylab("Peak location (%)") +
  xlab(NULL) +
  scale_fill_manual(
    values = txfeat_cols,
    guide = guide_legend(title = "Feature", reverse = TRUE)
  ) +
  ggtitle("Overlap of consensus R-loops with tx features")
plt


## Density plot of peak scores
consInt3 <- consInt2 %>%
  group_by(group) %>% 
  mutate(
    pct_cons = case_when(
      grepl(group, pattern = "dRNH") ~ 100 * ((score / 10) / ndrnh),
      grepl(group, pattern = "S9.6") ~ 100 * ((score / 10) / ns96)
    ),
    pct_rank = minmax(rank(pct_cons)),
  ) %>% 
  filter(group != "Shared")
pltdat <- consInt3 %>% distinct(name, .keep_all = TRUE) 
mu <- pltdat %>% group_by(group) %>% summarise(mu=mean(pct_cons))
plt <- ggplot(pltdat, aes(x = pct_cons, color = group, fill = group)) +
  geom_density(alpha=.3, adjust=2) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic(base_size = 18) +
  ylab("Frequency") + xlab("Conserved (%)") +
  guides(fill=guide_legend(title=NULL), color=guide_legend(title=NULL))
plt


## TSS / TES enrichment metaplot

# library(ChIPseeker)
# pks <- pltdat %>% 
#   filter(! group %in% c("S9.6-only", "S9.6-shared")) %>% 
#   group_by(group) %>% 
#   {
#     setNames(
#       group_split(.),
#       nm = group_keys(.)[[1]]
#     )
#   } %>% 
#   as.list() %>% 
#   lapply(GenomicRanges::makeGRangesFromDataFrame, keep.extra.columns = TRUE)
# 
# ## TSS
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# pms <- lapply(
#   pks, 
#   function(x) {
#     getTagMatrix(
#       peak = x, TxDb = txdb, 
#       upstream = 5000, downstream = 5000, 
#       type = "start_site", by = "gene", 
#       weightCol = "pct_cons",
#     )
#   }
# )
# 
# plotAvgProf(pms, xlim=c(-5000, 5000), free_y = F,  facet = "row", origin_label = "TSS")
# # tagHeatmap(pms, xlim=c(-5000, 5000), color=NULL)
# 
# ## TTS
# tms <- lapply(
#   pks, 
#   function(x) {
#     getTagMatrix(
#       peak = x, TxDb = txdb, 
#       upstream = 5000, downstream = 5000, 
#       type = "end_site", by = "gene", 
#       weightCol = "pct_cons",
#     )
#   }
# )
# 
# plotAvgProf(tms, xlim=c(-5000, 5000), free_y = F, facet = "row", origin_label = "TTS")
# # tagHeatmap(tms, xlim=c(-5000, 5000), color=NULL)
# 
### Genomic coverage plot
peakAnnoList <- lapply(pks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=T)

### Top 1000 pathway enrichment
geneLst <- consInt3 %>% 
  filter(! is.na(SYMBOL)) %>% 
  slice_max(order_by = pct_cons, n = 1000) %>% 
  {setNames(group_split(.), group_keys(.)[[1]])} %>%  
  lapply(pull, var = SYMBOL)
resRmd <- lapply(names(geneLst), function(groupNow) {
  message(groupNow)
  genesNow <- geneLst[[groupNow]]
  response <- httr::POST(  # Post request to enrichr based on https://maayanlab.cloud/Enrichr/help#api&q=1
    url = 'https://maayanlab.cloud/Enrichr/addList', 
    body = list(
      'list' = paste0(genesNow, collapse = "\n"),
      'description' = groupNow
    )
  )
  response <- jsonlite::fromJSON(httr::content(response, as = "text"))  # Collect response
  permalink <- paste0("https://maayanlab.cloud/Enrichr/enrich?dataset=",  # Create permalink
                      response$shortId[1])
  permalink
})
resRmd

## Make pathway enrichment plot
library(enrichR)
rllstenr <- pbapply::pblapply(
  geneLst,
  function(x) {
    genlst <- unique(unlist(strsplit(x, ",")))
    enrichR::enrichr(genes = genlst, databases = c("GO_Biological_Process_2021", "ChEA_2016",
                                                   "KEGG_2021_Human", "MSigDB_Hallmark_2020"))
  }
)

num_sel <- 6
db <- "KEGG_2021_Human"
db <- "ChEA_2016"
db <- "MSigDB_Hallmark_2020"
# db <- "GO_Biological_Process_2021"
# db <- "GO_Biological_Process_2021"
terms <- lapply(names(rllstenr), function(x) {
  rllstenr[[x]][[db]] %>%
    as_tibble() %>%
    # mutate(Combined.Score = as.numeric(scale(Combined.Score, center = FALSE))) %>%
    slice_max(Combined.Score, n = num_sel) %>%
    filter(row_number() <= num_sel) %>% pull(Term)
})
terms <- unique(unlist(terms))
plttbl <- lapply(names(rllstenr), function(x) {
  rllstenr[[x]][[db]] %>%
    as_tibble() %>%
    filter(Term %in% terms) %>%
    mutate(quant = x) %>%
    # mutate(Combined.Score = as.numeric(scale(Combined.Score, center = FALSE))) %>%
    select(Term, quant, combined_score=Combined.Score, padj=Adjusted.P.value)
}) %>% bind_rows()
barplt <- plttbl %>%
  mutate(
    combined_score = case_when(
      is.na(combined_score) | combined_score < 1 ~ 0,
      combined_score > 100 ~ 100,
      TRUE ~ combined_score
    )
  ) %>%
  arrange(quant, combined_score) %>% 
  group_by(quant) %>%
  mutate(
    Term = gsub(Term, pattern = " [0-9]+ .+", replacement = ""),
    Term = gsub(Term, pattern = " \\(GO.+", replacement = ""),
    Term = stringr::str_wrap(Term, width = 45),
  ) 
plt <- barplt %>%
  mutate(
    Term = factor(Term,
                  levels = unique(barplt$Term[barplt$combined_score != 0]))
  ) %>%
  filter(! quant %in% c("S9.6-only", "S9.6-shared")) %>% 
  filter(! is.na(Term)) %>%
  ggplot(aes(x = Term, color = -log10(padj), size=combined_score, y = quant)) +
  geom_point() +
  # geom_col(position = position_dodge(.2)) +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_viridis_c(direction = -1, option = "A", end = .9,
                        guide = guide_colorbar(title = "P adj (-log10)")) +
  ggtitle("TF target enrichment in consensus R-loop regions")
plt


### Chromatin States
states <- annFull[sapply(names(annFull), grepl, pattern="Encode_CREs")] %>% 
  bind_rows(.id = "group")

state_pks <- valr::bed_intersect(
  states, consInt3
)

statlst <- states %>% 
  group_by(group) %>% 
  {setNames(group_split(.), nm = group_keys(.)[[1]])}

pltdat <- consInt3 %>% 
  {setNames(group_split(.), nm = group_keys(.)[[1]])} %>% 
  lapply(function(x) {
    lapply(
      statlst, function(y) {
        valr::bed_fisher(x, y = y, genome = genome) %>% 
          mutate(
            state = y$group[1]
          )
      }
    ) %>% 
      bind_rows() %>% 
      mutate(
        group = x$group[1]
      )
  }) %>%
  bind_rows() %>% 
  mutate(
    p.value = ifelse(p.value == 0, .Machine$double.xmin, p.value),
    p.value = p.adjust(p.value),
    state = gsub(state, pattern = ".+__(.+)", replacement = "\\1"),
    state = case_when(
      state == "CTCF" ~ "CTCF-only",
      state == "enhD" ~ "Distal Enhancer",
      state == "enhP" ~ "Enhancer-Promoter",
      state == "prom" ~ "Promoter",
      state == "K4m3" ~ "H3K4me3"
    )
  ) %>% 
  pivot_wider(id_cols = group, names_from = state, values_from = estimate) %>% 
  column_to_rownames("group")

pltdat %>%
  log2() %>%
  t() %>% 
  pheatmap::pheatmap(
    color = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(200),
    angle_col = 45,
    main = "Enrichment of peaks within CREs",
    cluster_cols = F
  )

print(valr::bed_fisher(states, y = consInt3, genome = genome))


#### Analysis of dRNH-only peaks and enhancers #####

### Step 0: Get all the enhancers
if (! file.exists("data/genehancer_processed.rds")) {
  ghraw <- readxl::read_excel("data/GeneHancer_version_4-4.xlsx")
  attr <- ghraw %>% 
    pull(attributes) %>%
    str_split(pattern = ";") %>% 
    parallel::mclapply(function(x) {
      ghid <- x[grepl(x, pattern = "geneha")]
      genes <- x[grepl(x, pattern = "_gene")] %>%
        str_replace(pattern = ".+=", replacement = "")
      gene_scores <- x[grepl(x, pattern = "score")] %>%
        str_replace(pattern = ".+=", replacement = "") %>% 
        as.numeric()
      tibble(
        ID = gsub(ghid, pattern = ".+=", replacement = ""),
        genes, gene_scores
      ) %>% 
        group_by(ID) %>%
        nest()
    }, mc.cores = 44)
  d <- attr
  d <- split(d, ceiling(seq_along(d)/2000))
  d2 <- parallel::mclapply(d, bind_rows, mc.cores = 44)
  attr2 <- bind_rows(d2)
  gh <- bind_cols(ghraw, attr2) %>% 
    select(-attributes, name=ID) %>% 
    relocate(chrom, start, end, name, score, strand)
  saveRDS(gh, file = "data/genehancer_processed.rds")
} 
gh <- readRDS("data/genehancer_processed.rds")
ghfull <- gh %>% 
  unnest(cols = data)

## Find the distal enhancers
ghgr <- gh %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
ghanno <- ghgr %>% 
  annotatePeak(TxDb=txdb, tssRegion=c(-3000, 3000), verbose=T)
dEnh <- ghgr[which(ghanno@detailGenomicAnnotation$Intergenic),]

### Step 1: Establish further the relationship to enhancers for dRNH

pks3 <- pks %>% lapply(function(x) {x$score <- x$pct_cons; x})

## Overlap of dRNH-only & S9.6-only & all ENH
olrm <- ChIPpeakAnno::findOverlapsOfPeaks(pks$`dRNH-only`, pks$`dRNH-shared`)
ole <- ChIPpeakAnno::findOverlapsOfPeaks(ghgr, pks$`dRNH-only`, pks$`dRNH-shared`, pks$S9.6)
ChIPpeakAnno::makeVennDiagram(ole, NameOfPeaks = c("Enh", "dRNH-only", "dRNH-shared", "S9.6"))

gr2bed2 <- function(x) {
  as.data.frame(x) %>% 
    as_tibble() %>% 
    dplyr::rename(chrom=seqnames) 
}

enhstats <- lapply(pks, function(x) {
  y <- gr2bed2(x)
  ints <- valr::bed_intersect(y, gr2bed2(ghgr))
  nms <- unique(y$name)
  tibble(
    "names" = nms
  ) %>%
    mutate(inenh = ifelse(names %in% ints$name.x, "Enhancer", "non-Enhancer")) %>% 
    group_by(inenh) %>% 
    tally() %>% 
    mutate(
      prop = n / sum(n)
    )
}) %>% bind_rows(.id = "group")

enhstats %>% 
  mutate(group = factor(group, levels = rev(unique(.$group))),
         inenh = factor(inenh, levels = rev(unique(.$inenh)))) %>% 
  ggplot(aes(x = group, y = prop, fill=inenh)) +
  geom_col() +
  coord_flip() +
  theme_gray(base_size = 14) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  ggtitle("Proportion of S9.6 and dRNH peaks found at any enhancers")

plt <- ChIPpeakAnno::binOverFeature(
  pks3$`dRNH-shared`, pks3$`dRNH-only`, pks3$S9.6,
  annotationData=ghgr,
  radius=5000, nbins=100, FUN=length, errFun = 0,
  xlab="Distance from Enh (bp)", ylab="count", 
  main=c("dRNH-shared", "dRNH-only", "S9.6")
)
plt %>%
  as.data.frame() %>% 
  rownames_to_column("pos") %>% 
  as_tibble() %>% 
  mutate(pos = as.numeric(pos)) %>% 
  dplyr::rename("dRNH-shared"=2, "dRNH-only"=3, "S9.6"=4) %>% 
  pivot_longer(cols = ! contains("pos")) %>% 
  ggplot(aes(x = pos, y = value, color = name)) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=.25) +
  geom_point() +
  geom_line() +
  ggtitle("Consensus peak pileup around all enhancers") +
  ylab("Peak density") +
  xlab("Distance to any enhancer (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("dRNH-only" = "#e2a3e3", "dRNH-shared" = "#ccb6cc", "S9.6" = "#82d0e8")) 

## Overlap of dRNH-only & S9.6-only & dENH
# olde <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, pks$`dRNH-only`, pks$`dRNH-shared`, pks$S9.6)
# ChIPpeakAnno::makeVennDiagram(olde, NameOfPeaks = c("Distal Enh", "dRNH-only", "dRNH-shared", "S9.6"))
enhstats <- lapply(pks, function(x) {
  y <- gr2bed2(x)
  ints <- valr::bed_intersect(y, gr2bed2(dEnh))
  nms <- unique(y$name)
  tibble(
    "names" = nms
  ) %>%
    mutate(inenh = ifelse(names %in% ints$name.x, "Enhancer", "non-Enhancer")) %>% 
    group_by(inenh) %>% 
    tally() %>% 
    mutate(
      prop = n / sum(n)
    )
}) %>% bind_rows(.id = "group")

enhstats %>% 
  mutate(group = factor(group, levels = rev(unique(.$group))),
         inenh = factor(inenh, levels = rev(unique(.$inenh)))) %>% 
  ggplot(aes(x = group, y = prop, fill=inenh)) +
  geom_col() +
  coord_flip() +
  theme_gray(base_size = 14) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  ggtitle("Proportion of S9.6 and dRNH peaks found at distal enhancers")
plt <- ChIPpeakAnno::binOverFeature(
  pks3$`dRNH-shared`, pks3$`dRNH-only`, pks3$S9.6,
  annotationData=dEnh,
  radius=5000, nbins=100, FUN=length, errFun = 0,
  xlab="Distance from dEnh (bp)", ylab="count", 
  main=c("dRNH-shared", "dRNH-only", "S9.6")
)
plt %>%
  as.data.frame() %>% 
  rownames_to_column("pos") %>% 
  as_tibble() %>% 
  mutate(pos = as.numeric(pos)) %>% 
  dplyr::rename("dRNH-shared"=2, "dRNH-only"=3, "S9.6"=4) %>% 
  pivot_longer(cols = ! contains("pos")) %>% 
  ggplot(aes(x = pos, y = value, color = name)) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=.25) +
  geom_point() +
  geom_line() +
  ggtitle("Consensus peak pileup around distal enhancers") +
  ylab("Peak density") +
  xlab("Distance to distal enhancer (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("dRNH-only" = "#e2a3e3", "dRNH-shared" = "#ccb6cc", "S9.6" = "#82d0e8")) 
  

## Overlap of intergenic dRNH-only & S9.6-only & dENH
pks4 <- lapply(seq(length(pks)), function(i) {
  pks[[i]][which(peakAnnoList[[i]]@detailGenomicAnnotation$Intergenic),]
})
names(pks4) <- names(pks)
olde2 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, pks4$`dRNH-only`, pks4$`dRNH-shared`, pks4$S9.6)
ChIPpeakAnno::makeVennDiagram(olde2, NameOfPeaks = c("Distal Enh", "dRNH-only (IG)", "dRNH-shared (IG)", "S9.6 (IG)"))
ChIPpeakAnno::binOverFeature(
  pks4$`dRNH-shared`, pks4$`dRNH-only`, pks4$S9.6,
  annotationData=dEnh,
  radius=5000, nbins=100, FUN=length, errFun = 0,
  xlab="Distance from dEnh (bp)", ylab="count", 
  main=c("dRNH-shared (IG)", "dRNH-only (IG)","S9.6 (IG)")
)


### Step 2: What are the genes which these dENH's interact with?
topdRNH <- pks4$`dRNH-only`
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, topdRNH)
ChIPpeakAnno::makeVennDiagram(olde3, NameOfPeaks = c("Distal-Enhancers", "dRNH-only"), margin=.1, cat.pos=c(270, 90))
dENH_dRNH <- olde3$overlappingPeaks$`dEnh///topdRNH`$name
ghfull %>% 
  filter(
    name %in% {{ dENH_dRNH }},
    gene_scores > 30
  ) -> dd
unique(dd$genes) -> genesNow
genesNow
groupNow <- "dRNH-bound dEnh targets"
resRmd <- lapply(names(geneLst), function(groupNow) {
  message(groupNow)
  genesNow <- geneLst[[groupNow]]
  response <- httr::POST(  # Post request to enrichr based on https://maayanlab.cloud/Enrichr/help#api&q=1
    url = 'https://maayanlab.cloud/Enrichr/addList', 
    body = list(
      'list' = paste0(genesNow, collapse = "\n"),
      'description' = groupNow
    )
  )
  response <- jsonlite::fromJSON(httr::content(response, as = "text"))  # Collect response
  permalink <- paste0("https://maayanlab.cloud/Enrichr/enrich?dataset=",  # Create permalink
                      response$shortId[1])
  permalink
})
resRmd

write_csv(tibble(genesNow), file = "analyses/extension_of_dRNH_S96_analysis/dEnh_dRNH-only_gene_targets.csv")

#### Step #3: Tissue-specific enhancer analysis ####

dat <- RLHub::rlbase_samples()

### CUTLL1

ct <- "CUTLL1"

## Load in the data
ctl <- dat %>%
  filter(
    ip_type == "dRNH",
    genome == "hg38",
    tissue == ct,
    label == "POS",
    prediction == "POS",
    numPeaks > 5000
  ) 
grs <- lapply(ctl$rlsample, function(x) {
  GenomicRanges::GRanges(RLSeq::RLRangesFromRLBase(x))
})
names(grs) <- ctl$rlsample

## Analyze feature distribution
pal <- lapply(
  grs, annotatePeak, TxDb=txdb,
  tssRegion=c(-3000, 3000), verbose=T
)

## Find overlaps with enhancers

# Get intergenic ranges
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)

# Test enrichment within dEnh
olde <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, gris[[1]], gris[[2]])
ChIPpeakAnno::makeVennDiagram(olde, NameOfPeaks = c("Enh", names(gris)[1], names(gris)[2]))
ChIPpeakAnno::binOverFeature(
  gris[[1]], gris[[2]],
  annotationData=dEnh,
  radius=5000, nbins=100, FUN=length, errFun = 0,
  xlab="Distance from Enh (bp)", ylab="count",
  main=c(names(gris)[1], names(gris)[2])
)

## Find the specific Enh
# Get the intergenic peaks from both iPSC samples & find overlap
gri <- GenomicRanges::reduce(do.call("c", unlist(gris, use.names = F)), with.revmap=T)
gri <- gri[which(sapply(gri$revmap, length) > 1),]
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, gri)
ChIPpeakAnno::makeVennDiagram(olde3)

# For these overlapping peaks, what are they?
dENH_gri <- olde3$overlappingPeaks$`dEnh///gri`
ghfull %>% 
  filter(
    name %in% dENH_gri$name,
    gene_scores > 10
  ) -> dd
unique(dd$genes) -> genesNow
genesNow
eres <- enrichr(genesNow, databases = "CellMarker_Augmented_2021")
num_sel <- 6
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>% 
  filter(row_number() <= num_sel) %>% pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms) %>%
  select(Term, combined_score=Combined.Score, padj=Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>% 
  mutate(
    Term = factor(Term,
                  levels = unique(Term))
  ) %>%
  filter(! is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("CUTLL1 dRNH-enhancer gene targets", subtitle = "Enrichment in CellMarker Database")
plt


### iPSCs
ct <- "iPSCs"

## Load in the data
ctl <- dat %>%
  filter(
    ip_type == "dRNH",
    genome == "hg38",
    tissue == ct,
    label == "POS",
    prediction == "POS",
    numPeaks > 5000
  ) 
grs <- lapply(ctl$rlsample, function(x) {
  GenomicRanges::GRanges(RLSeq::RLRangesFromRLBase(x))
})
names(grs) <- ctl$rlsample

## Analyze feature distribution
pal <- lapply(
  grs, annotatePeak, TxDb=txdb,
  tssRegion=c(-3000, 3000), verbose=T
)

## Find overlaps with enhancers

# Get intergenic ranges
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)

## Find the specific Enh
# Get the intergenic peaks from both iPSC samples & find overlap
gri <- GenomicRanges::reduce(do.call("c", unlist(gris, use.names = F)), with.revmap=T)
gri <- gri[which(sapply(gri$revmap, length) > 1),]
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, gri)
ChIPpeakAnno::makeVennDiagram(olde3)

# For these overlapping peaks, what are they?
dENH_gri <- olde3$overlappingPeaks$`dEnh///gri`
ghfull %>% 
  filter(
    name %in% dENH_gri$name,
    gene_scores > 20
  ) -> dd
unique(dd$genes) -> genesNow
genesNow
eres <- enrichr(genesNow, databases = "CellMarker_Augmented_2021")
num_sel <- 6
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>% 
  filter(row_number() <= num_sel) %>% pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms) %>%
  select(Term, combined_score=Combined.Score, padj=Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>% 
  mutate(
    Term = factor(Term,
                  levels = unique(Term))
  ) %>%
  filter(! is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer gene targets", subtitle = "Enrichment in CellMarker Database")
plt


#### Step #4: Chromatin State Analysis (poised vs active enhancers) ####

fls <- list(
  "HEK293" = "https://www.encodeproject.org/files/ENCFF476TTU/@@download/ENCFF476TTU.bed.gz",
  "HCT116" = "https://www.encodeproject.org/files/ENCFF513PJK/@@download/ENCFF513PJK.bed.gz",
  "iPSCs" = "https://www.encodeproject.org/files/ENCFF115RIR/@@download/ENCFF115RIR.bed.gz",
  "iPSCs2" = "https://www.encodeproject.org/files/ENCFF519CTX/@@download/ENCFF519CTX.bed.gz",
  "iPSCs3" = "https://www.encodeproject.org/files/ENCFF676VUR/@@download/ENCFF676VUR.bed.gz"
)

pks <- lapply(
  fls, read_tsv, 
  col_names = c(
    "chrom", "start", "end", "name", 
    "score", "strand", "start2", "end2", "color"
  )
)
pks <- lapply(
  pks, function(x) {
    filter(x, str_detect(name, "Enh")) %>% 
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  }
)

# Look at overlap between dENH and predicted enhancers
dd <- GenomicRanges::reduce(do.call("c", unlist(pks, use.names = FALSE)))
ol1 <- ChIPpeakAnno::findOverlapsOfPeaks(dd, dEnh)
ChIPpeakAnno::makeVennDiagram(ol1)


### iPSCs

ct <- "iPSCs"

## Load in the data
ctl <- dat %>%
  filter(
    ip_type == "dRNH",
    genome == "hg38",
    tissue == ct,
    label == "POS",
    prediction == "POS",
    numPeaks > 5000
  ) 
grs <- lapply(ctl$rlsample, function(x) {
  gr <- GenomicRanges::GRanges(RLSeq::RLRangesFromRLBase(x))
  # gr <- gr[gr$qval > 3,]
  gr
})
names(grs) <- ctl$rlsample

## Analyze feature distribution
pal <- lapply(
  grs, annotatePeak, TxDb=txdb,
  tssRegion=c(-3000, 3000), verbose=T
)

## Find overlaps with enhancers

# Get intergenic ranges
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)
gri <- GenomicRanges::reduce(do.call("c", unlist(gris, use.names = F)), with.revmap=T)
gri <- gri[which(sapply(gri$revmap, length) > 1),]

## Find the specific Enh
# Get the intergenic peaks from both iPSC samples & find overlap
pki <- pks$iPSCs %>% 
  annotatePeak(TxDb = txdb, tssRegion = c(-3000, 3000), verbose = TRUE)
pki <- pks$iPSCs[pki@detailGenomicAnnotation$Intergenic,]
pki <- pki[! pki$name %in% c("EnhG1", "EnhG2"), ]

# Transfer labels to the dEnh from GH
olpe <- ChIPpeakAnno::findOverlapsOfPeaks(
  pki, dEnh[dEnh$score > .5,]
)
ChIPpeakAnno::makeVennDiagram(olpe)
colnames(olpe$overlappingPeaks$`pki///dEnh`)[7] <- "pName"
dEnh2 <- olpe$overlappingPeaks$`pki///dEnh` %>% 
  select(name, pName) %>% 
  as_tibble() %>% 
  group_by(name) %>% 
  mutate(n=n()) %>% 
  filter(n==1) %>% 
  select(-n) %>% 
  inner_join(as.data.frame(dEnh)) %>% 
  as.data.frame() %>% 
  ChIPpeakAnno::toGRanges()


# Analyze dEnh and Peaks
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh2, gri)
ChIPpeakAnno::makeVennDiagram(olde3)
d1 <- olde3$overlappingPeaks$`dEnh2///gri`$pName %>% 
  tibble() %>% 
  dplyr::rename(enh=1) %>% 
  group_by(enh) %>% 
  tally() %>% 
  mutate(`dRNH-accessible` = n / sum(n)) %>% 
  select(-n)

olde3$all.peaks$dEnh2$pName %>%
  tibble() %>% 
  dplyr::rename(enh=1) %>% 
  group_by(enh) %>% 
  tally() %>% 
  mutate(`Total dEnh` = n / sum(n)) %>% 
  select(-n) %>% 
  inner_join(d1) %>% 
  pivot_longer(cols = ! matches("^enh")) %>% 
  ggplot(aes(fill = enh, y = value, x = name)) +
  geom_col() +
  coord_flip() +
  xlab(NULL) +
  ylab("Proportion of distal enhancers") +
  ggtitle(
    paste0(ct, " distal enhancer profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_gray(base_size = 14)


# For these overlapping peaks, what are they?
dENH_gri <- olde3$overlappingPeaks$`dEnh2///gri`
ghfull %>% 
  filter(
    name %in% dENH_gri$peaks1[dENH_gri$pName %in% c("EnhBiv")],
    gene_scores > 10
  ) -> dd
unique(dd$genes) -> genesNow
genesNow
groupNow <- paste0(ct, " - Enh analysis")
response <- httr::POST(  # Post request to enrichr based on https://maayanlab.cloud/Enrichr/help#api&q=1
  url = 'https://maayanlab.cloud/Enrichr/addList', 
  body = list(
    'list' = paste0(genesNow, collapse = "\n"),
    'description' = groupNow
  )
)
response <- jsonlite::fromJSON(httr::content(response, as = "text"))  # Collect response
permalink <- paste0("https://maayanlab.cloud/Enrichr/enrich?dataset=",  # Create permalink
                    response$shortId[1])
permalink


eres <- enrichr(genesNow, databases = "CellMarker_Augmented_2021")
num_sel <- 6
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>% 
  filter(row_number() <= num_sel) %>% pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms) %>%
  select(Term, combined_score=Combined.Score, padj=Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>% 
  mutate(
    Term = factor(Term,
                  levels = unique(Term))
  ) %>%
  filter(! is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in CellMarker Database")
plt

eres <- enrichr(genesNow, databases = "ChEA_2016")
num_sel <- 6
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>% 
  filter(row_number() <= num_sel) %>% pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms) %>%
  select(Term, combined_score=Combined.Score, padj=Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>% 
  mutate(
    Term = factor(Term,
                  levels = unique(Term))
  ) %>%
  filter(! is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer gene targets", subtitle = "Enrichment in ChEA Database")
plt


eres <- enrichr(genesNow, databases = "ARCHS4_TFs_Coexp")
num_sel <- 6
eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel)
terms <- eres[[1]] %>%
  as_tibble() %>%
  slice_min(P.value, n = num_sel) %>% 
  filter(row_number() <= num_sel) %>% pull(Term)
terms <- unique(unlist(terms))
terms
plttbl <- eres[[1]] %>%
  as_tibble() %>%
  filter(Term %in% terms, P.value < 0.01) %>%
  select(Term, combined_score=Combined.Score, padj=Adjusted.P.value)
plt <- plttbl %>%
  arrange(combined_score) %>% 
  mutate(
    Term = factor(Term,
                  levels = unique(Term))
  ) %>%
  filter(! is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in ARCHS4 TF Database")
plt

#### Step #5: Analyze relationship with eRNA 

chain <-  'data/hg38_to_hg19.chain.gz'
if (! file.exists(gsub(chain, pattern="\\.gz", replacement=""))) {
  download.file("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz", destfile = chain)
  R.utils::gunzip(chain)
}

### iPSCs

gro1 <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3271nnn/GSM3271001/suppl/GSM3271001_iPSC_11104_A.bedGraph.gz",
  col_names = c("chrom", "start", "end", "score")
)
gro2 <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3271nnn/GSM3271002/suppl/GSM3271002_iPSC_11104_B.bedGraph.gz",
  col_names = c("chrom", "start", "end", "score")
)

int <- valr::bed_intersect(gro1, gro2)
gro <- int %>% 
  mutate(score = score.x + score.y) %>% 
  select(chrom, start=start.x, end=end.y, score)


chn <- rtracklayer::import.chain(gsub(chain, pattern = "\\.gz", replacement = ""))
dEnh2hg19 <- rtracklayer::liftOver(dEnh2, chain = chn) %>% unlist()
dEnh2hg19 <- dEnh2hg19[! duplicated(names(dEnh2hg19)),]

dEnh2_min <- as.data.frame(dEnh2hg19) %>% 
  as_tibble() %>% 
  mutate(name = names(dEnh2hg19)) %>% 
  select(chrom=seqnames, start, end, name)
int <- valr::bed_intersect(gro, dEnh2_min)
intGRO <- int %>%
  group_by(name.y) %>% 
  summarise(GRO = sum(score.x)) %>% 
  dplyr::rename(name=name.y)
# Get RPKM
dEnh_min2 <- dEnh2_min %>% mutate(width = end - start)
intGRO <- inner_join(intGRO, select(dEnh_min2, name, width)) %>% 
  mutate(
    RPKM = GRO /
      (
        (width/1000) *
          (sum(GRO)/1000000)
      ),
    RPKM = log2(RPKM + 1)
 ) %>% select(-width)
dEnh3 <- as.data.frame(dEnh2hg19) %>% 
  as_tibble() %>% 
  mutate(name = names(dEnh2hg19),
         name2=name) %>% 
  inner_join(intGRO) %>% 
  as.data.frame() %>% 
  ChIPpeakAnno::toGRanges() %>% 
  unique()

# hg19 the GRI
grihg19 <- rtracklayer::liftOver(gri, chain = chn) %>% unlist()

# Overlap with iPSC peaks
olde4 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh3, grihg19)
ChIPpeakAnno::makeVennDiagram(olde4)
d1 <- olde4$overlappingPeaks$`dEnh3///grihg19`$RPKM %>%  
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "dRNH-accessible"
  )

d2 <- olde4$all.peaks$dEnh3$RPKM %>% 
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "Total dEnh"
  )


bind_rows(d1, d2) %>% 
  filter(value > 0) %>%
  ggplot(aes(y = value, x = cond)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(
    method.args = list(alternative="greater"),
    comparisons = list(c("dRNH-accessible", "Total dEnh")),
    label = "p.format"
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_gray(base_size = 14)



d1 <- olde4$overlappingPeaks$`dEnh3///grihg19`$RPKM[olde4$overlappingPeaks$`dEnh3///grihg19`$pName == "EnhBiv"] %>% 
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "dRNH-accessible"
  )

d2 <- olde4$all.peaks$dEnh3$RPKM[ olde4$all.peaks$dEnh3$pName == "EnhBiv"] %>%
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "Total dEnh"
  )


bind_rows(d1, d2) %>% 
  filter(value > 0) %>%
  ggplot(aes(y = value, x = cond)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(
    method.args = list(alternative="greater"),
    comparisons = list(c("dRNH-accessible", "Total dEnh"))
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile (Bivalent)"),
    subtitle = "dRNH-accessible vs total bivalent enhancer population"
  ) +
  theme_gray(base_size = 14)


### CUTLL1
ct <- "CUTLL1"
gro1 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_DMSO_forward.bw",
  set_strand = "*"
)
gro2 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_DMSO_reverse.bw",
  set_strand = "*"
) %>% mutate(score=-1*score)
gro3 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_00h_forward.bw",
  set_strand = "*"
)
gro4 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_00h_reverse.bw",
  set_strand = "*"
) %>% mutate(score=-1*score)

gro <- valr::bed_merge(bind_rows(list(
  gro1, 
  gro2, 
  gro3, 
  gro4
)), score = sum(score))

int <- valr::bed_intersect(gro, dEnh2_min)
intGRO <- int %>%
  group_by(name.y) %>% 
  summarise(GRO = sum(score.x)) %>% 
  dplyr::rename(name=name.y)
# Get RPKM
intGRO
intGRO <- inner_join(intGRO, select(dEnh_min2, name, width)) %>% 
  mutate(
    RPKM = GRO /
      (
        (width/1000) *
          (sum(GRO)/1000000)
      ),
    RPKM = log2(RPKM + 1)
  ) %>% select(-width)
dEnh3 <- as.data.frame(dEnh2hg19) %>% 
  as_tibble() %>% 
  mutate(name = names(dEnh2hg19),
         name2=name) %>% 
  inner_join(intGRO) %>% 
  as.data.frame() %>% 
  ChIPpeakAnno::toGRanges() %>% 
  unique()

# hg19 the GRI
grihg19 <- rtracklayer::liftOver(gri, chain = chn) %>% unlist()

# Overlap with iPSC peaks
olde4 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh3, grihg19)
ChIPpeakAnno::makeVennDiagram(olde4)
d1 <- olde4$overlappingPeaks$`dEnh3///grihg19`$RPKM %>% #[olde4$overlappingPeaks$`dEnh3///gri`$pName == "EnhBiv"] %>% 
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "dRNH-accessible"
  )

d2 <- olde4$all.peaks$dEnh3$RPKM %>% #[ olde4$all.peaks$dEnh3$pName == "EnhBiv"] %>%
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "Total dEnh"
  )


bind_rows(d1, d2) %>%
  ggplot(aes(y = value, x = cond)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(
    method.args = list(alternative="greater"),
    comparisons = list(c("dRNH-accessible", "Total dEnh")),
    label = "p.format"
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_gray(base_size = 14)



d1 <- olde4$overlappingPeaks$`dEnh3///grihg19`$RPKM[olde4$overlappingPeaks$`dEnh3///grihg19`$pName == "EnhBiv"] %>% 
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "dRNH-accessible"
  )

d2 <- olde4$all.peaks$dEnh3$RPKM[ olde4$all.peaks$dEnh3$pName == "EnhBiv"] %>%
  tibble() %>% 
  dplyr::rename(value=1) %>% 
  mutate(
    cond = "Total dEnh"
  )


bind_rows(d1, d2) %>% 
  filter(value > 0) %>%
  ggplot(aes(y = value, x = cond)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(
    method.args = list(alternative="greater"),
    comparisons = list(c("dRNH-accessible", "Total dEnh"))
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile (Bivalent)"),
    subtitle = "dRNH-accessible vs total bivalent enhancer population"
  ) +
  theme_gray(base_size = 14)






