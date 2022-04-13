## This is an extension of the analysis from the original paper
## Purpose is to characterize S96 and dRNH peaks and dig deeper into them

library(tidyverse)
library(ChIPseeker)
library(enrichR)

#### MAGIC ####

ndrnh <- 56 # Number of dRNH passing QC filters
ns96 <- 191 # Number of S9.6 passing QC filters
num_sel <- 6
annFull <- RLHub::annots_full_hg38()
resloc <- "results/New_Figres"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

###############

dir.create(resloc, showWarnings = FALSE)

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
) %>%
  dplyr::rename(chrom = SEQNAME, start = GENESEQSTART, end = GENESEQEND) %>%
  filter(!grepl(GENEID, pattern = "LRG.+")) %>%
  select(-GENEID) %>%
  mutate(chrom = paste0("chr", chrom)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

## Get the dRNH-only and S96-only along with shared
int <- valr::bed_intersect(s96_cons, drnh_cons)
cons2 <- list(
  "dRNH-only" = drnh_cons[!drnh_cons$name %in% int$name.y, ],
  "S9.6-only" = s96_cons[!s96_cons$name %in% int$name.x, ],
  "dRNH-shared" = drnh_cons[drnh_cons$name %in% int$name.y, ],
  "S9.6-shared" = s96_cons[s96_cons$name %in% int$name.x, ],
  "S9.6" = s96_cons
)


## Density plot of peak sizes
minmax <- function(x) (x - min(x)) / (max(x) - min(x))
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
mu <- pltdat %>%
  group_by(group2) %>%
  summarise(mu = mean(width))
mu
plt <- ggplot(pltdat, aes(y = width, x = group2, fill = group2)) +
  geom_violin(width = .6) +
  geom_boxplot(width = .2) +
  scale_y_log10() +
  theme_classic(base_size = 18) +
  ylab("Peak width (log bp)") +
  xlab(NULL) +
  ggtitle("Consensus peak width distribution") +
  guides(
    fill = guide_legend(title = NULL), 
    color = guide_legend(title = NULL)
  ) +
  scale_fill_manual(
    values = c(
      "dRNH-only" = "#e2a3e3",
      "dRNH-shared" = "#af91bf", 
      "S9.6" = "#82d0e8"
    )
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("dRNH-only", "dRNH-shared"))
  ) +
  ggpubr::rremove("legend")
plt
ggsave(plt, filename = file.path(resloc, "consensus_peak_width.png"))

#### Step #2: Proportion test
## What prop of dRNH-detected R-loops are also detected by S9.6?
## What prop of S9.6-detected R-loops are also detected by dRNH?
## Check with pie charts or Bar charts
plt <- pltdat %>%
  filter(group != "S9.6") %>%
  group_by(group) %>%
  tally() %>%
  mutate(
    group2 = gsub(group, pattern = "\\-.+", replacement = ""),
    group3 = gsub(group, pattern = ".+\\-", replacement = ""),
    group3 = ifelse(group3 == "only", "Unshared", "Shared"),
    group3 = factor(group3, levels = (unique(group3)))
  ) %>%
  group_by(group2) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = group2, y = prop, fill = group3)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  ggtitle("RLRegions shared vs unshared")
plt
ggsave(plt, filename = file.path(resloc, "prop_consensus_shared_barchart.png"))


#### Step #3: Analysis of Tx Features
### What are the features overlapping consensus peaks in each group?
### How does this relate to peak size/strength?
cons3 <- cons2[!names(cons2) %in% c("S9.6-only", "S9.6-shared")]
## Compare samples around genomic features
consInt2 <- lapply(cons3, function(x) {
  # Call summits
  x <- x %>%
    dplyr::rename(
      so = start,
      eo = end
    ) %>%
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
                           ifelse("Exon" %in% type.y, "Exon", "Intron")
                    )
             )
      )
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
    n = 100 * n / sum(n),
    type = factor(type, levels = rev(c(
      "TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "Intergenic"
    )))
  ) %>%
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
ggsave(plt, filename = file.path(resloc, "consensus_peak_tx_overlap.png"))

## Wrangle peak lists
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
pks <- pltdat %>%
  filter(!group %in% c("S9.6-only", "S9.6-shared")) %>%
  group_by(group) %>%
  {
    setNames(
      group_split(.),
      nm = group_keys(.)[[1]]
    )
  } %>%
  as.list() %>%
  lapply(GenomicRanges::makeGRangesFromDataFrame, keep.extra.columns = TRUE)


### Top 1000 pathway enrichment
geneLst <- consInt3 %>%
  filter(!is.na(SYMBOL)) %>%
  slice_max(order_by = pct_cons, n = 1000) %>%
  {
    setNames(group_split(.), group_keys(.)[[1]])
  } %>%
  lapply(pull, var = SYMBOL)

## Make pathway enrichment plot
rllstenr <- pbapply::pblapply(
  geneLst,
  function(x) {
    genlst <- unique(unlist(strsplit(x, ",")))
    enrichR::enrichr(
      genes = genlst, databases = c(
        "GO_Biological_Process_2021", "ChEA_2016",
        "KEGG_2021_Human", "MSigDB_Hallmark_2020"
      )
    )
  }
)

db <- "KEGG_2021_Human"
db <- "ChEA_2016"
db <- "MSigDB_Hallmark_2020"
# db <- "GO_Biological_Process_2021"
# db <- "GO_Biological_Process_2021"
terms <- lapply(names(rllstenr), function(x) {
  rllstenr[[x]][[db]] %>%
    as_tibble() %>%
    slice_max(Combined.Score, n = num_sel) %>%
    filter(row_number() <= num_sel) %>%
    pull(Term)
})
terms <- unique(unlist(terms))
plttbl <- lapply(names(rllstenr), function(x) {
  rllstenr[[x]][[db]] %>%
    as_tibble() %>%
    filter(Term %in% terms) %>%
    mutate(quant = x) %>%
    select(
      Term,
      quant, 
      combined_score = Combined.Score, 
      padj = Adjusted.P.value
    )
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
    Term = factor(
      Term,
      levels = unique(barplt$Term[barplt$combined_score != 0])
    )
  ) %>%
  filter(!quant %in% c("S9.6-only", "S9.6-shared")) %>%
  filter(!is.na(Term)) %>%
  ggplot(
    aes(x = Term, color = -log10(padj), size = combined_score, y = quant)
  ) +
  geom_point() +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_viridis_c(
    direction = -1, option = "A", end = .9,
    guide = guide_colorbar(title = "P adj (-log10)")
  ) +
  ggtitle("Consensus region pathway enrichment")
plt
ggsave(plt, filename = file.path(resloc, "consensus_pathway_enrichment.png"))

### Chromatin States
states <- annFull[sapply(names(annFull), grepl, pattern = "Encode_CREs")] %>%
  bind_rows(.id = "group")

state_pks <- valr::bed_intersect(
  states, consInt3
)

statlst <- states %>%
  group_by(group) %>%
  {
    setNames(group_split(.), nm = group_keys(.)[[1]])
  }

pltdat <- consInt3 %>%
  {
    setNames(group_split(.), nm = group_keys(.)[[1]])
  } %>%
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

plt <- pltdat %>%
  log2() %>%
  t() %>%
  pheatmap::pheatmap(
    color = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(200),
    angle_col = 45,
    main = "Enrichment of peaks within CREs",
    cluster_cols = F
  )
png(file.path(resloc, "CCRE_heatmap.png"))
plt
dev.off()


#### Analysis of dRNH-only peaks and enhancers #####

### Step 0: Get all the enhancers
if (!file.exists("data/genehancer_processed.rds")) {
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
  d <- split(d, ceiling(seq_along(d) / 2000))
  d2 <- parallel::mclapply(d, bind_rows, mc.cores = 44)
  attr2 <- bind_rows(d2)
  gh <- bind_cols(ghraw, attr2) %>%
    select(-attributes, name = ID) %>%
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
  annotatePeak(TxDb = txdb, tssRegion = c(-3000, 3000), verbose = T)
dEnh <- ghgr[which(ghanno@detailGenomicAnnotation$Intergenic), ]

### Step 1: Establish further the relationship to enhancers for dRNH
pks3 <- pks %>% lapply(function(x) {
  x$score <- x$pct_cons
  x
})

## Overlap of dRNH-only & S9.6-only & all ENH
gr2bed2 <- function(x) {
  as.data.frame(x) %>%
    as_tibble() %>%
    dplyr::rename(chrom = seqnames)
}

enhstats <- lapply(pks, function(x) {
  y <- gr2bed2(x)
  ints <- valr::bed_intersect(y, gr2bed2(ghgr))
  nms <- unique(y$name)
  tibble(
    "names" = nms
  ) %>%
    mutate(
      inenh = ifelse(names %in% ints$name.x, "Enhancer", "non-Enhancer")
    ) %>%
    group_by(inenh) %>%
    tally() %>%
    mutate(
      prop = n / sum(n)
    )
}) %>% bind_rows(.id = "group")

plt <- enhstats %>%
  mutate(
    group = factor(group, levels = rev(unique(.$group))),
    inenh = factor(inenh, levels = rev(unique(.$inenh)))
  ) %>%
  ggplot(aes(x = group, y = prop, fill = inenh)) +
  geom_col() +
  coord_flip() +
  theme_gray(base_size = 14) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  ggtitle("Consensus peaks found at enhancers")
plt
ggsave(plt, filename = file.path(resloc, "consensus_pks_all_enhancers.png"))

pltd <- ChIPpeakAnno::binOverFeature(
  pks3$`dRNH-shared`, pks3$`dRNH-only`, pks3$S9.6,
  annotationData = ghgr,
  radius = 5000, nbins = 100, FUN = length, errFun = 0,
  xlab = "Distance from Enh (bp)", ylab = "count",
  main = c("dRNH-shared", "dRNH-only", "S9.6")
)
plt <- pltd %>%
  as.data.frame() %>%
  rownames_to_column("pos") %>%
  as_tibble() %>%
  mutate(pos = as.numeric(pos)) %>%
  dplyr::rename("dRNH-shared" = 2, "dRNH-only" = 3, "S9.6" = 4) %>%
  pivot_longer(cols = !contains("pos")) %>%
  ggplot(aes(x = pos, y = value, color = name)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .25) +
  geom_point() +
  geom_line() +
  ggtitle("Consensus peak pileup around all enhancers") +
  ylab("Peak density") +
  xlab("Distance to any enhancer (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_color_manual(
    values = c(
      "dRNH-only" = "#e2a3e3", 
      "dRNH-shared" = "#ccb6cc", 
      "S9.6" = "#82d0e8"
    )
  )
plt
ggsave(plt, filename = file.path(resloc, "consensus_pileup_all_enh.png"))

## Overlap of dRNH-only & S9.6-only & dENH
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

plt <- enhstats %>%
  mutate(
    group = factor(group, levels = rev(unique(.$group))),
    inenh = gsub(inenh, pattern = "Enh", replacement = "distal-enh"),
  ) %>%
  mutate(
    inenh = factor(inenh, levels = rev(unique(.$inenh)))
  ) %>% 
  ggplot(aes(x = group, y = prop, fill = inenh)) +
  geom_col() +
  coord_flip() +
  theme_gray(base_size = 14) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  ggtitle("Proportion of S9.6 and dRNH peaks found at distal enhancers")
plt
ggsave(plt, filename = file.path(resloc, "consensus_pks_denh.png"))

pltd <- ChIPpeakAnno::binOverFeature(
  pks3$`dRNH-shared`, pks3$`dRNH-only`, pks3$S9.6,
  annotationData = dEnh,
  radius = 5000, nbins = 100, FUN = length, errFun = 0,
  xlab = "Distance from dEnh (bp)", ylab = "count",
  main = c("dRNH-shared", "dRNH-only", "S9.6")
)
plt <- pltd %>%
  as.data.frame() %>%
  rownames_to_column("pos") %>%
  as_tibble() %>%
  mutate(pos = as.numeric(pos)) %>%
  dplyr::rename("dRNH-shared" = 2, "dRNH-only" = 3, "S9.6" = 4) %>%
  pivot_longer(cols = !contains("pos")) %>%
  ggplot(aes(x = pos, y = value, color = name)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .25) +
  geom_point() +
  geom_line() +
  ggtitle("Consensus peak pileup around distal enhancers") +
  ylab("Peak density") +
  xlab("Distance to distal enhancer (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_color_manual(
    values = c(
      "dRNH-only" = "#e2a3e3", 
      "dRNH-shared" = "#ccb6cc", 
      "S9.6" = "#82d0e8"
    )
  )
plt
ggsave(plt, filename = file.path(resloc, "consensus_pileup_dEnh.png"))

## Overlap of intergenic dRNH-only & S9.6-only & dENH
peakAnnoList <- lapply(pks, annotatePeak,
                       TxDb = txdb,
                       tssRegion = c(-3000, 3000), verbose = T
)
pks4 <- lapply(seq(length(pks)), function(i) {
  pks[[i]][which(peakAnnoList[[i]]@detailGenomicAnnotation$Intergenic), ]
})
names(pks4) <- names(pks)

topdRNH <- pks4$`dRNH-only`
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, topdRNH)
png(file.path(resloc, "Venn_dRNHonlt_dEnh.png"))
ChIPpeakAnno::makeVennDiagram(
  olde3, 
  NameOfPeaks = c("Distal-Enhancers", "Intergenic dRNH-only"),
  margin = .1, 
  fill = c("grey", "#e2a3e3"),
  cat.pos = c(240, 120)
)
dev.off()

#### Step #2: Tissue-specific enhancer analysis ####

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
  grs, annotatePeak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000), verbose = T
)

## Find overlaps with enhancers

# Get intergenic ranges
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)

# Test enrichment within dEnh
pltd <- ChIPpeakAnno::binOverFeature(
  gris[[1]], gris[[2]], gris[[3]], gris[[4]], gris[[5]], gris[[6]],
  annotationData = dEnh,
  radius = 5000, nbins = 100, FUN = length, errFun = 0,
  xlab = "Distance from Enh (bp)", ylab = "count",
  main = c(names(gris)[1], names(gris)[2])
)
plt <- pltd %>%
  as.data.frame() %>%
  rownames_to_column("pos") %>%
  as_tibble() %>%
  mutate(pos = as.numeric(pos)) %>%
  dplyr::rename(
    "SRX10484271" = 2, 
    "SRX10484272" = 3, 
    "SRX10484273" = 4,
    "SRX10484274" = 5,
    "SRX10484275" = 6,
    "SRX10484276" = 7
  ) %>%
  pivot_longer(cols = !contains("pos")) %>%
  ggplot(aes(x = pos, y = value, color = name)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .25) +
  geom_point() +
  geom_line() +
  ggtitle("Peak pileup around distal enhancers", subtitle = "CUTLL1 MapR (dRNH) peaks") +
  ylab("Peak density") +
  xlab("Distance to distal enhancer (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2")
plt
ggsave(plt, filename = file.path(resloc, "CUTLL1_pileup_dEnh.png"))

## Find the specific Enh
# Get the intergenic peaks from both iPSC samples & find overlap
gri <- GenomicRanges::reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, gri)
# PCT overlapping
length(unique(olde3$overlappingPeaks$`dEnh///gri`$peaks2)) / length(unique(olde3$all.peaks$gri))
ChIPpeakAnno::makeVennDiagram(
  olde3,
  NameOfPeaks = c("Distal Enhancers", "Intergenic dRNH sites"),
  fill=c("gray", "#0BB490"), margin=.1
)

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
  geom_col(fill="#0BB490") +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("CUTLL1 dRNH-enhancer gene targets", subtitle = "Enrichment in CellMarker Database")
plt
ggsave(plt, filename = file.path(resloc, "Cellmarker_CUTLL1_dENH.png"))


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
  grs, annotatePeak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000), verbose = T
)

## Find overlaps with enhancers

# Get intergenic ranges
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)

# Test enrichment within dEnh
pltd <- ChIPpeakAnno::binOverFeature(
  gris[[1]], gris[[2]], gris[[3]], gris[[4]],
  annotationData = dEnh,
  radius = 5000, nbins = 100, FUN = length, errFun = 0,
  xlab = "Distance from Enh (bp)", ylab = "count"
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
  ggtitle("Peak pileup around distal enhancers", subtitle = "iPSCs MapR (dRNH) peaks") +
  ylab("Peak density") +
  xlab("Distance to distal enhancer (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2")
plt
ggsave(plt, filename = file.path(resloc, "iPSCs_pileup_dEnh.png"))

## Find the specific Enh
# Get the intergenic peaks from both iPSC samples & find overlap
gri <- reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]
olde3 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh, gri)
# PCT olap
length(unique(olde3$overlappingPeaks$`dEnh///gri`$peaks2)) / length(unique(olde3$all.peaks$gri))
ChIPpeakAnno::makeVennDiagram(
  olde3,
  NameOfPeaks = c("Distal Enhancers", "Intergenic dRNH sites"),
  fill=c("gray", "#F58776"), margin=.1
)

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
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer gene targets", subtitle = "Enrichment in CellMarker Database")
plt


#### Step #3: Chromatin State Analysis (poised vs active enhancers) ####

fls <- list(
  "HEK293" = "https://www.encodeproject.org/files/ENCFF476TTU/@@download/ENCFF476TTU.bed.gz",
  "HCT116" = "https://www.encodeproject.org/files/ENCFF513PJK/@@download/ENCFF513PJK.bed.gz",
  "iPSCs" = "https://www.encodeproject.org/files/ENCFF115RIR/@@download/ENCFF115RIR.bed.gz"
)

# Get track for hg19 chromHMM
# iPSChg19 <- "https://www.encodeproject.org/files/ENCFF197OWG/@@download/ENCFF197OWG.bed.gz"
# infile <- "tmp/iPSC_chromhmm_enhbiv.bed"
# infile2 <- "tmp/iPSC_chromhmm_enhbiv.sort.bed"
# outfile <- "tmp/iPSC_chromhmm_enhbiv.sort.bb"
# iPSChg19 <- valr::read_bed12(iPSChg19)
# iPSChg19 %>% filter(name == "EnhBiv") %>% 
#   write_tsv("tmp/iPSC_chromhmm_enhbiv.bed", col_names = F)
# 
# system(paste0("tail -n +2 ", infile, " | sort -k1,1 -k2,2n > ", infile2))
# system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedToBigBed ", infile2, " tmp/chrom.sizes ", outfile))
# aws.s3::put_object(file = outfile, object = "misc/iPSC_hg19_enhbig.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)

# TODO: Add in the DRIPc data + DRIP data
# TODO: Look at the other places where you have this effect (look for biv enhancer pileup)
# TODO: To get the U pattern, do a signal metaplot around bivalent enhancer clusters -- should should GRO-SEq U shape



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

# Look at overlap between dENH and predicted enhancers
dd <- reduce(do.call("c", unlist(pkschrom, use.names = FALSE)))
ol1 <- ChIPpeakAnno::findOverlapsOfPeaks(dd, dEnh)

### iPSCs

ct <- "iPSCs"

## Find the specific Enh
# Get the intergenic peaks from both iPSC samples & find overlap
pki <- pkschrom$iPSCs %>%
  annotatePeak(TxDb = txdb, tssRegion = c(-3000, 3000), verbose = TRUE)
pki <- pkschrom$iPSCs[pki@detailGenomicAnnotation$Intergenic, ]
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
  geom_col() +
  coord_flip() +
  xlab(NULL) +
  ylab("Proportion of distal enhancers") +
  ggtitle(
    paste0(ct, " distal enhancer profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_gray(base_size = 14)
plt
ggsave(plt, filename = file.path(resloc, "iPSC_dENH_types_prop.png"))

# For these overlapping peaks, what are they?
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
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in CellMarker Database")
plt
ggsave(plt, filename = file.path(resloc, "iPSC_dRNH-dEnh_Biv_cellmarker.png"))

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
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in ChEA Database")
plt
ggsave(plt, filename = file.path(resloc, "iPSC_dRNH-dEnh_Biv_chea.png"))

eres <- enrichr(genesNow, databases = "ARCHS4_TFs_Coexp")
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
    Term = factor(Term,
                  levels = unique(Term)
    )
  ) %>%
  filter(!is.na(Term)) %>%
  ggplot(aes(x = Term, y = combined_score)) +
  geom_col(fill="#F58776") +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer (Biv) gene targets", subtitle = "Enrichment in ARCHS4 TF Database")
plt
ggsave(plt, filename = file.path(resloc, "iPSC_dRNH-dEnh_Biv_archs4.png"))


#### Analyze bivalent enhancer clusters

ct <- "iPSCs"
chain <- "data/hg38_to_hg19.chain.gz"
chn <- rtracklayer::import.chain(gsub(chain, pattern = "\\.gz", replacement = ""))

# Calculate clusters for all enhancers
pkbiv <- pkschrom$iPSCs 
pkbivhg19 <- rtracklayer::liftOver(pkbiv, chain = chn) %>% unlist()
pkbivhg19 <- unique(pkbivhg19)
bigbiv <- (pkbivhg19 + 5000) %>% GenomicRanges::reduce(with.revmap=T)
len <- sapply(bigbiv$revmap, length)
bigbiv[len>1,] %>% 
  rtracklayer::export(
    con = "tmp/enh_clusters.bed"
  )

system(
  paste0(
    '~/miniconda3/condabin/conda run -n rlbaseData computeMatrix scale-regions ',
    '-S "analyses/extension_of_dRNH_S96_analysis/iPSC_GROSeq_comb.bw" ',
    '"analyses/extension_of_dRNH_S96_analysis/iPSC_MapR_comb.bw" ',
    'tmp/GSE136857_iPSC_IDR.bw ',
    # 'tmp/ENCFF645JYS.bigWig tmp/ENCFF402HIK.bigWig ',
    'tmp/rlregions_dRNH.hg19.bw.bw tmp/rlregions_S96.hg19.bw.bw -R ',
    'tmp/bivenh_clusters.bed tmp/enh_clusters.bed -m 15000 -a 15000 -b 15000 -bs 300 -o ',
    'tmp/enh_clusters.mat.gz -p 44'
  )
)

n <- 5

mat <- read_tsv("tmp/enh_clusters.mat.gz", skip = 1, col_names = FALSE)
meta <- read_lines("tmp/enh_clusters.mat.gz")[1]
meta <- gsub(meta, pattern = "@", replacement = "")
meta <- jsonlite::parse_json(meta)

colnames(mat)[1:6] <- c("chrom", "start", "end", "name", "score", "strand")
colnames(mat)[7:(6+(n*150))] <- paste0(
  rep(
    c(
      "GRO", 
      "MapR",
      "ATAC",
      # "H3K27me3",
      # "H3K4me1",
      "dRNH",
      "S96"
    ), each=150
  ),
  "__",
  rep(1:50, n*3),
  "__",
  rep(1:150, n),
  "__",
  rep(rep(c("a", "m", "b"), each=50), n)
)

mlong <- mat %>% pivot_longer(cols = contains("__"), names_to = "NAMES")
lpat <- "(.+)__(.+)__(.+)__(.+)"
mlong <- mlong %>% 
  mutate(
    ip=gsub(NAMES, pattern = lpat, replacement = "\\1"),
    pos = gsub(NAMES, pattern = lpat, replacement = "\\2"),
    posful =  gsub(NAMES, pattern = lpat, replacement = "\\3"),
    posgroup = gsub(NAMES, pattern = lpat, replacement = "\\4"),
    value = as.numeric(case_when(value == "NaN" ~ 0, TRUE ~ value))
  )

sums <- mlong %>%
  group_by(ip) %>% 
  summarise(val=sum(value))

mlong2 <- inner_join(mlong, sums)
mlong3 <- mlong2 %>%
  mutate(value_scale = round(1000000*(value / val), digits = 6))
pltdat <- mlong3 %>% 
  group_by(ip, posful) %>% 
  summarise(value_scale2 = mean(value_scale))
pltdat %>% 
  ggplot(aes(x = posful, y = value_scale2, color=ip, group=ip)) +
  geom_line()

mwide <- mlong3 %>% 
  select(-value, -ip, -pos, -posful, -posgroup, -val) %>% 
  pivot_wider(names_from = NAMES, values_from = value_scale)

write_tsv(mwide, file = "tmp/enh_clusters2.mat", col_names = FALSE)
system("zcat tmp/enh_clusters.mat.gz | head -n 1 > tmp/enh_clusters3.mat && cat tmp/enh_clusters2.mat >> tmp/enh_clusters3.mat && gzip -f tmp/enh_clusters3.mat")
system(
  paste0(
    '~/miniconda3/condabin/conda run -n rlbaseData plotHeatmap --matrixFile ',
    'tmp/enh_clusters3.mat.gz --outFileName tmp/enh_clusters.png ',
    '-z "Bivalent enhancers" "All Enhancers" --samplesLabel "iPSC GRO-Seq" ',
    '"iPSC MapR" "iPSC ATAC-Seq" "dRNH consensus" "S9.6 consensus" ',
    '--startLabel	"start" --endLabel "end" --heatmapWidth 6 --colorMap "plasma"'
  )
)

# Calculate clusters for bivalent enhancers
pkchrom <- pkschrom$iPSCs 
pkbiv <- pkchrom[pkchrom$name == "EnhBiv",]
pkbivhg19 <- rtracklayer::liftOver(pkbiv, chain = chn) %>% unlist()
pkbivhg19 <- unique(pkbivhg19)
bigbiv <- (pkbivhg19 + 5000) %>% GenomicRanges::reduce(with.revmap=T)
len <- sapply(bigbiv$revmap, length)
bigbiv[len>1,] %>% 
  rtracklayer::export(
    con = "tmp/bivenh_clusters.bed"
  )

system(
  paste0(
    '~/miniconda3/condabin/conda run -n rlbaseData computeMatrix scale-regions ',
    '-S "analyses/extension_of_dRNH_S96_analysis/iPSC_GROSeq_comb.bw" ',
    '"analyses/extension_of_dRNH_S96_analysis/iPSC_MapR_comb.bw" ',
    'tmp/GSE136857_iPSC_IDR.bw tmp/ENCFF645JYS.bigWig tmp/ENCFF402HIK.bigWig ',
    'tmp/rlregions_dRNH.hg19.bw.bw tmp/rlregions_S96.hg19.bw.bw -R ',
    '"tmp/bivenh_clusters.bed" -m 15000 -a 15000 -b 15000 -bs 300 -o ',
    'tmp/bivenh_clusters.mat.gz -p 44'
  )
)

# TODO: Add DRIP into 
# TODO: Add 10KB bumpers


# Convert S9.6 and dRNH consensus to hg19




# Transfer labels to the dEnh from GH
olpe <- ChIPpeakAnno::findOverlapsOfPeaks(
  pki, dEnh[dEnh$score > .5, ]
)



##### Step #4: Analyze relationship with eRNA

chain <- "data/hg38_to_hg19.chain.gz"
chn <- rtracklayer::import.chain(gsub(chain, pattern = "\\.gz", replacement = ""))
if (!file.exists(gsub(chain, pattern = "\\.gz", replacement = ""))) {
  download.file(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz", 
    destfile = chain
  )
  R.utils::gunzip(chain)
}

#### iPSCs

ct <- "iPSCs"

# Read in GRO-Seq
gro1 <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3271nnn/GSM3271001/suppl/GSM3271001_iPSC_11104_A.bedGraph.gz",
  col_names = c("chrom", "start", "end", "score")
)
gro2 <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3271nnn/GSM3271002/suppl/GSM3271002_iPSC_11104_B.bedGraph.gz",
  col_names = c("chrom", "start", "end", "score")
)

# Get intersect
int <- valr::bed_intersect(gro1, gro2)
gro <- int %>%
  mutate(score = score.x + score.y) %>%
  select(chrom, start = start.x, end = end.y, score)

# Convert dRNH data to hg19
dEnh2hg19 <- rtracklayer::liftOver(dEnh2, chain = chn) %>% unlist()
dEnh2hg19 <- dEnh2hg19[!duplicated(names(dEnh2hg19)), ]

# Get minimal form of hg19 dENH data
dEnh2_min <- as.data.frame(dEnh2hg19) %>%
  as_tibble() %>%
  mutate(name = names(dEnh2hg19)) %>%
  select(chrom = seqnames, start, end, name)

# Intersect with GRO-Seq to find eRNA signal at dENH
int <- valr::bed_intersect(gro, dEnh2_min)
intGRO <- int %>%
  group_by(name.y) %>%
  summarise(GRO = sum(score.x)) %>%
  dplyr::rename(name = name.y)

# Get eRNA RPKM
dEnh_min2 <- dEnh2_min %>% mutate(width = end - start)
intGRO <- inner_join(intGRO, select(dEnh_min2, name, width)) %>%
  mutate(
    RPKM = GRO /
      (
        (width / 1000) *
          (sum(GRO) / 1000000)
      ),
    RPKM = log2(RPKM + 1)
  ) %>%
  select(-width)
dEnh3 <- as.data.frame(dEnh2hg19) %>%
  as_tibble() %>%
  mutate(
    name = names(dEnh2hg19),
    name2 = name
  ) %>%
  inner_join(intGRO) %>%
  as.data.frame() %>%
  ChIPpeakAnno::toGRanges() %>%
  unique()

# Convert the iPSC MapR intergenic data to hg19
grihg19 <- rtracklayer::liftOver(gri, chain = chn) %>% unlist()

# Overlap RPKM eRNA dataset with MapR data
olde4 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh3, grihg19)
ChIPpeakAnno::makeVennDiagram(olde4)

# Calculate eRNA levels at MapR-accessible dENH vs all dENH
d1 <- olde4$overlappingPeaks$`dEnh3///grihg19`$RPKM %>%
  tibble() %>%
  dplyr::rename(value = 1) %>%
  mutate(
    cond = "dRNH-accessible dEnh"
  )
d2 <- olde4$all.peaks$dEnh3$RPKM %>%
  tibble() %>%
  dplyr::rename(value = 1) %>%
  mutate(
    cond = "Total dEnh"
  )
plt <- bind_rows(d1, d2) %>%
  filter(value > 0) %>%
  ggplot(aes(y = value, x = cond, fill=cond)) +
  geom_boxplot(width=.6) +
  ggpubr::stat_compare_means(
    method.args = list(alternative = "greater"),
    comparisons = list(c("dRNH-accessible dEnh", "Total dEnh")),
    label = "p.format"
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(
    values = c(
      "dRNH-accessible dEnh" = "#F58776",
      "Total dEnh" = "grey"
    )
  ) +
  ggpubr::rremove("legend")

plt
ggsave(plt, filename = file.path(resloc, "iPSC_eRNA_expression_dENH.png"))

# ### Prep genome browser session for showing this
# cs19 <- "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"
# download.file(cs19, destfile = "tmp/hg19.chrom.sizes")
# gen <- valr::read_genome("tmp/hg19.chrom.sizes")
# gen %>% mutate(start = 0, end = size) -> dd
# wins <- valr::bed_makewindows(x = dd, win_size = 50)
# gm1 <- valr::bed_map(wins, gro1, score2 = mean(score))
# gm2 <- valr::bed_map(wins, gro2, score2 = mean(score))
# GRO <- select(gm1, chrom, start, end)
# GRO$score <- (gm1$score2 + gm2$score2) / 2
# infile <- "analyses/extension_of_dRNH_S96_analysis/iPSC_GROSeq_comb.bdg"
# infile2 <- "analyses/extension_of_dRNH_S96_analysis/iPSC_GROSeq_comb.sort.bdg"
# outfile <- "analyses/extension_of_dRNH_S96_analysis/iPSC_GROSeq_comb.bw"
# GRO %>% write_tsv(infile, col_names = FALSE)
# system(paste0("tail -n +2 ", infile, " | sort -k1,1 -k2,2n > ", infile2))
# system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedGraphToBigWig ", infile2, " tmp/chrom.sizes ", outfile))
# aws.s3::put_object(file = outfile, object = "misc/iPSC_GROSeq_comb.bw", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
# 
# 
# ## Do the iPSC MapR
# bws <- lapply(
#   ctl$rlsample, function(x) {
#     message(x)
#     bw <- paste0("../RLBase-data/rlbase-data/rlpipes-out/coverage/", x, "_hg38.bw")
#     bw <- valr::read_bigwig(bw, set_strand = "*")
#     bw
#   }
# )
# 
# # Lift over
# bws2 <- lapply(seq(bws), function(i) {
#   message(i)
#   x <- bws[[i]]
#   gr <- x %>%
#     GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
#   rtracklayer::liftOver(gr, chain = chn) %>%
#     unlist() %>%
#     as_tibble() %>%
#     rename(chrom = seqnames)
# })
# 
# # Map to windows
# bws3 <- parallel::mclapply(seq(bws2), function(i) {
#   message(i)
#   x <- bws2[[i]] %>%
#     select(-width) %>%
#     relocate(score, .before = strand)
#   valr::bed_map(wins, x, score2 = mean(score))
# }, mc.cores = 4)
# 
# # Combine
# infile <- "analyses/extension_of_dRNH_S96_analysis/iPSC_MapR_comb.bdg"
# infile2 <- "analyses/extension_of_dRNH_S96_analysis/iPSC_MapR_comb.sort.bdg"
# outfile <- "analyses/extension_of_dRNH_S96_analysis/iPSC_MapR_comb.bw"
# ipsrl <- bws3[[1]]
# ipsrl$score <- (bws3[[1]]$score2 + bws3[[2]]$score2 + bws3[[3]]$score2) / 3
# ipsrl %>%
#   filter(!is.na(score)) %>%
#   select(chrom, start, end, score) %>%
#   write_tsv(col_names = F)
# system(paste0("tail -n +2 ", infile, " | sort -k1,1 -k2,2n > ", infile2))
# system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedGraphToBigWig ", infile2, " tmp/chrom.sizes ", outfile))
# aws.s3::put_object(file = outfile, object = "misc/iPSC_MapR_comb.bw", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
# 
# 
# ## Upload the dRNH-only big bed
# download.file("https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes", destfile = "tmp/chrom.sizes")
# drnhhg19 <- cons2$`dRNH-only` %>%
#   mutate(score = score / 10) %>%
#   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
#   rtracklayer::liftOver(chain = chn) %>%
#   unlist() %>%
#   as_tibble() %>%
#   select(-width) %>%
#   relocate(name, score, .before = strand) %>%
#   mutate(strand = ".") %>%
#   unique()
# drnhhg19 %>% write_tsv("analyses/extension_of_dRNH_S96_analysis/dRNHonly.narrowPeak", col_names = F)
# infile <- "analyses/extension_of_dRNH_S96_analysis/dRNHonly.narrowPeak"
# infile2 <- "analyses/extension_of_dRNH_S96_analysis/dRNHonly.sort.narrowPeak"
# outfile <- "analyses/extension_of_dRNH_S96_analysis/dRNHonly.bb"
# 
# 
# # Re-upload hg19 versions of all consensus peaks
# sapply(c("S96", "dRNH"), function(x) {
#   download.file("https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes", destfile = "tmp/hg19.chrom.sizes")
#   infile <- paste0("../RLBase-data/rlbase-data/rlregions/rlregions_", x, ".narrowPeak")
#   infile2 <- paste0("analyses/extension_of_dRNH_S96_analysis/rlregions_", x, ".hg19.narrowPeak")
#   dat <- valr::read_narrowpeak(infile)
#   dat <- dat %>%
#     mutate(score = score / 10) %>%
#     GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
#     rtracklayer::liftOver(chain = chn) %>%
#     unlist() %>%
#     as_tibble() %>%
#     select(-width) %>%
#     relocate(name, score, .before = strand) %>%
#     mutate(strand = ".") %>%
#     unique()
#   write_tsv(dat, infile2, col_names = FALSE)
#   infile3 <- paste0("analyses/extension_of_dRNH_S96_analysis/rlregions_", x, ".hg19.sort.narrowPeak")
#   outfile <- paste0("analyses/extension_of_dRNH_S96_analysis/rlregions_", x, ".hg19.bb")
#   system(paste0("tail -n +2 ", infile2, " | sort -k1,1 -k2,2n > ", infile3))
#   system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedToBigBed -type=bed6+4 -as=tmp/bigNarrowPeak.as ", infile3, " tmp/hg19.chrom.sizes ", outfile))
# })
# aws.s3::put_object(file = "analyses/extension_of_dRNH_S96_analysis/rlregions_S96.hg19.bb", object = "misc/rlregions_S96.hg19.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
# aws.s3::put_object(file = "analyses/extension_of_dRNH_S96_analysis/rlregions_dRNH.hg19.bb", object = "misc/rlregions_dRNH.hg19.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
# 
# 
# 
# drnhhg19 <- cons2$`dRNH-only` %>%
#   mutate(score = score / 10) %>%
#   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
#   rtracklayer::liftOver(chain = chn) %>%
#   unlist() %>%
#   as_tibble() %>%
#   select(-width) %>%
#   relocate(name, score, .before = strand) %>%
#   mutate(strand = ".") %>%
#   unique()
# drnhhg19 %>% write_tsv("analyses/extension_of_dRNH_S96_analysis/dRNHonly.narrowPeak", col_names = F)
# infile <- "analyses/extension_of_dRNH_S96_analysis/dRNHonly.narrowPeak"
# infile2 <- "analyses/extension_of_dRNH_S96_analysis/dRNHonly.sort.narrowPeak"
# outfile <- "analyses/extension_of_dRNH_S96_analysis/dRNHonly.bb"
# 
# 
# 
# download.file("https://genome.ucsc.edu/goldenPath/help/examples/bigNarrowPeak.as", destfile = "tmp/bigNarrowPeak.as")
# system(paste0("tail -n +2 ", infile, " | sort -k1,1 -k2,2n > ", infile2))
# system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedToBigBed -type=bed6+4 -as=tmp/bigNarrowPeak.as ", infile2, " tmp/chrom.sizes ", outfile))
# aws.s3::put_object(file = outfile, object = "misc/dRNHonly.bb", bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
# 
# 
# 
# 


#### CUTLL1
ct <- "CUTLL1"
gro1 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_DMSO_forward.bw",
  set_strand = "*"
)
gro2 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_DMSO_reverse.bw",
  set_strand = "*"
) %>% mutate(score = -1 * score)
gro3 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_00h_forward.bw",
  set_strand = "*"
)
gro4 <- valr::read_bigwig(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115894/suppl/GSE115894_00h_reverse.bw",
  set_strand = "*"
) %>% mutate(score = -1 * score)

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
  dplyr::rename(name = name.y)

# Get RPKM
intGRO <- inner_join(intGRO, select(dEnh_min2, name, width)) %>%
  mutate(
    RPKM = GRO /
      (
        (width / 1000) *
          (sum(GRO) / 1000000)
      ),
    RPKM = log2(RPKM + 1)
  ) %>%
  select(-width)
dEnh3 <- as.data.frame(dEnh2hg19) %>%
  as_tibble() %>%
  mutate(
    name = names(dEnh2hg19),
    name2 = name
  ) %>%
  inner_join(intGRO) %>%
  as.data.frame() %>%
  ChIPpeakAnno::toGRanges() %>%
  unique()

# Convert the MapR intergenic peaks to hg19
grihg19 <- rtracklayer::liftOver(gri, chain = chn) %>% unlist()

# Overlap with iPSC peaks
olde4 <- ChIPpeakAnno::findOverlapsOfPeaks(dEnh3, grihg19)
ChIPpeakAnno::makeVennDiagram(olde4)
d1 <- olde4$overlappingPeaks$`dEnh3///grihg19`$RPKM %>% 
  tibble() %>%
  dplyr::rename(value = 1) %>%
  mutate(
    cond = "dRNH-bound dEnh"
  )
d2 <- olde4$all.peaks$dEnh3$RPKM %>% 
  tibble() %>%
  dplyr::rename(value = 1) %>%
  mutate(
    cond = "Total dEnh"
  )
plt <- bind_rows(d1, d2) %>%
  ggplot(aes(y = value, x = cond, fill=cond)) +
  geom_boxplot(width=.6) +
  ggpubr::stat_compare_means(
    method.args = list(alternative = "greater"),
    comparisons = list(c("dRNH-bound dEnh", "Total dEnh")),
    label = "p.format"
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(
    values = c(
      "dRNH-bound dEnh" = "#0BB490",
      "Total dEnh" = "gray"
    )
  ) +
  ggpubr::rremove("legend")
plt
ggsave(plt, filename = file.path(resloc, "CUTLL1_eRNA_expression_dENH.png"))

##### ANALYSIS of eCLiP/ChIP datasets #####

rlrf <- RLHub::feat_enrich_rlregions()

annbind <- annFull[grepl(names(annFull), pattern = "eCLiP|ChIP")]

pks2 <- lapply(cons2[c("dRNH-only", "S9.6-only")], function(x) {
  # Call summits
  x <- x %>%
    dplyr::rename(
      so = start,
      eo = end
    ) %>%
    mutate(
      start = so + peak - 250,
      end = so + peak + 250
    )
  xint <- x %>%
    relocate(end, .after = chrom) %>% 
    relocate(start, .after = chrom)
  return(xint)
})

### Full
reslst1 <- lapply(pks2, function(x) {
  rlr <- RLSeq::RLRanges(
    peaks = GenomicRanges::makeGRangesFromDataFrame(x),
    genome = "hg38", mode = "DRIP"
  )
  rlr <- RLSeq::featureEnrich(annotations = annbind, object = rlr)
  rlr@metadata$results@featureEnrichment
})
dd <- bind_rows(reslst1, .id = "group") %>%
  group_by(group) %>%
  mutate(
    anno = paste0(gsub(db, pattern = "RBP_", replacement = ""), "__", type),
    stat_rl = -log10(pval_fisher_rl) * stat_fisher_rl,
    stat_fisher_rl = log2(stat_fisher_rl + 1),
    stat_shuf = -log10(pval_fisher_shuf) * stat_fisher_shuf,
    enrichment = stat_rl - stat_shuf,
    enrichment = case_when(enrichment < 0 ~ 0, TRUE ~ enrichment),
    enrichment = log2(enrichment + 1),
    enrichment = case_when(enrichment > 20 ~ 20, TRUE ~ enrichment),
    ip = gsub(db, pattern = "RBP_", replacement = "")
  ) %>%
  select(group, ip, anno, enrichment = enrichment) %>%
  pivot_wider(id_cols = c(anno, ip), names_from = group, values_from = enrichment) %>%
  mutate(anno = gsub(anno, pattern = ".+__", replacement = ""))
dd %>% 
  ggplot(aes(x = `S9.6-only`, y = `dRNH-only`, color = ip, label = anno)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
  geom_point() +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0, 15)) +
  theme_gray(base_size = 14) +
  ggtitle("All")
# scale_x_continuous(limits = c(0, 7)) +
# scale_y_continuous(limits = c(0, 7)) 
# ggrepel::geom_label_repel()



### Intron
reslst2 <- lapply(pks2, function(x) {
  pk <- GenomicRanges::makeGRangesFromDataFrame(x)
  ann <- annotatePeak(
    pk, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = TRUE
  )
  rlr <- RLSeq::RLRanges(
    peaks = pk[ann@detailGenomicAnnotation$Intron,],
    genome = "hg38", mode = "DRIP"
  )
  rlr <- RLSeq::featureEnrich(annotations = annbind, object = rlr)
  rlr@metadata$results@featureEnrichment
})
bind_rows(reslst2, .id = "group") %>%
  group_by(group) %>%
  mutate(
    anno = paste0(gsub(db, pattern = "RBP_", replacement = ""), "__", type),
    stat_rl = -log10(pval_fisher_rl) * stat_fisher_rl,
    stat_fisher_rl = log2(stat_fisher_rl + 1),
    stat_shuf = -log10(pval_fisher_shuf) * stat_fisher_shuf,
    enrichment = stat_rl - stat_shuf,
    enrichment = case_when(enrichment < 0 ~ 0, TRUE ~ enrichment),
    enrichment = log2(enrichment + 1),
    enrichment = case_when(enrichment > 20 ~ 20, TRUE ~ enrichment),
    ip = gsub(db, pattern = "RBP_", replacement = "")
  ) %>%
  select(group, ip, anno, enrichment = enrichment) %>%
  pivot_wider(id_cols = c(anno, ip), names_from = group, values_from = enrichment) %>%
  mutate(anno = gsub(anno, pattern = ".+__", replacement = "")) %>% 
  ggplot(aes(x = `S9.6-only`, y = `dRNH-only`, color = ip, label = anno)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
  geom_point() +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0, 15)) +
  theme_gray(base_size = 14) +
  ggtitle("Intron")
# scale_x_continuous(limits = c(0, 7)) +
# scale_y_continuous(limits = c(0, 7)) 
# ggrepel::geom_label_repel()



### Promoter
reslst3 <- lapply(pks2, function(x) {
  pk <- GenomicRanges::makeGRangesFromDataFrame(x)
  ann <- annotatePeak(
    pk, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = TRUE
  )
  rlr <- RLSeq::RLRanges(
    peaks = pk[ann@detailGenomicAnnotation$Promoter,],
    genome = "hg38", mode = "DRIP"
  )
  rlr <- RLSeq::featureEnrich(annotations = annbind, object = rlr)
  rlr@metadata$results@featureEnrichment
})
bind_rows(reslst3, .id = "group") %>%
  group_by(group) %>%
  mutate(
    anno = paste0(gsub(db, pattern = "RBP_", replacement = ""), "__", type),
    stat_rl = -log10(pval_fisher_rl) * stat_fisher_rl,
    stat_fisher_rl = log2(stat_fisher_rl + 1),
    stat_shuf = -log10(pval_fisher_shuf) * stat_fisher_shuf,
    enrichment = stat_rl - stat_shuf,
    enrichment = case_when(enrichment < 0 ~ 0, TRUE ~ enrichment),
    enrichment = log2(enrichment + 1),
    enrichment = case_when(enrichment > 20 ~ 20, TRUE ~ enrichment),
    ip = gsub(db, pattern = "RBP_", replacement = "")
  ) %>%
  select(group, ip, anno, enrichment = enrichment) %>%
  pivot_wider(id_cols = c(anno, ip), names_from = group, values_from = enrichment) %>%
  mutate(anno = gsub(anno, pattern = ".+__", replacement = "")) %>% 
  ggplot(aes(x = `S9.6-only`, y = `dRNH-only`, color = ip, label = anno)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
  geom_point() +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0, 15)) +
  theme_gray(base_size = 14) +
  ggtitle("Promoter")




### Distal Intergenic
reslst4 <- lapply(pks2, function(x) {
  pk <- GenomicRanges::makeGRangesFromDataFrame(x)
  ann <- annotatePeak(
    pk, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = TRUE
  )
  rlr <- RLSeq::RLRanges(
    peaks = pk[ann@detailGenomicAnnotation$distal_intergenic,],
    genome = "hg38", mode = "DRIP"
  )
  rlr <- RLSeq::featureEnrich(annotations = annbind, object = rlr)
  rlr@metadata$results@featureEnrichment
})
dd2 <- bind_rows(reslst4, .id = "group") %>%
  group_by(group) %>%
  mutate(
    anno = paste0(gsub(db, pattern = "RBP_", replacement = ""), "__", type),
    stat_rl = -log10(pval_fisher_rl) * stat_fisher_rl,
    stat_fisher_rl = log2(stat_fisher_rl + 1),
    stat_shuf = -log10(pval_fisher_shuf) * stat_fisher_shuf,
    enrichment = stat_rl - stat_shuf,
    enrichment = case_when(enrichment < 0 ~ 0, TRUE ~ enrichment),
    enrichment = log2(enrichment + 1),
    enrichment = case_when(enrichment > 20 ~ 20, TRUE ~ enrichment),
    ip = gsub(db, pattern = "RBP_", replacement = "")
  ) %>%
  select(group, ip, anno, enrichment = stat_fisher_rl) %>%
  pivot_wider(id_cols = c(anno, ip), names_from = group, values_from = enrichment) %>%
  mutate(anno = gsub(anno, pattern = ".+__", replacement = ""))
dd2 %>% 
  ggplot(aes(x = `S9.6-only`, y = `dRNH-only`, color = ip, label = anno)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
  geom_point() +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0, 15)) +
  theme_gray(base_size = 14) +
  ggtitle("Intergenic")



### Genic
reslst5 <- lapply(pks2, function(x) {
  pk <- GenomicRanges::makeGRangesFromDataFrame(x)
  ann <- annotatePeak(
    pk, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = TRUE
  )
  rlr <- RLSeq::RLRanges(
    peaks = pk[ann@detailGenomicAnnotation$genic,],
    genome = "hg38", mode = "DRIP"
  )
  rlr <- RLSeq::featureEnrich(annotations = annbind, object = rlr)
  rlr@metadata$results@featureEnrichment
})
dd3 <- bind_rows(reslst5, .id = "group") %>%
  group_by(group) %>%
  mutate(
    anno = paste0(gsub(db, pattern = "RBP_", replacement = ""), "__", type),
    stat_rl = -log10(pval_fisher_rl) * stat_fisher_rl,
    stat_fisher_rl = log2(stat_fisher_rl + 1),
    stat_shuf = -log10(pval_fisher_shuf) * stat_fisher_shuf,
    enrichment = stat_rl - stat_shuf,
    enrichment = case_when(enrichment < 0 ~ 0, TRUE ~ enrichment),
    enrichment = log2(enrichment + 1),
    enrichment = case_when(enrichment > 20 ~ 20, TRUE ~ enrichment),
    ip = gsub(db, pattern = "RBP_", replacement = "")
  ) %>%
  select(group, ip, anno, enrichment = stat_fisher_rl) %>%
  pivot_wider(id_cols = c(anno, ip), names_from = group, values_from = enrichment) %>%
  mutate(anno = gsub(anno, pattern = ".+__", replacement = ""))
dd3 %>% 
  ggplot(aes(x = `S9.6-only`, y = `dRNH-only`, color = ip, label = anno)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
  geom_point() +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(limits = c(0, 15)) +
  theme_gray(base_size = 14) +
  ggtitle("Intergenic")


### Comp
dd3 %>% 
  pivot_longer(cols = contains("-")) %>% 
  ggplot(aes(x = name, y = value, color = ip)) +
  geom_jitter(width = .15) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("dRNH-only", "S9.6-only")
    ), method = "t.test"
  ) +
  facet_grid(cols = vars(ip)) +
  theme_gray(13) +
  xlab(NULL) +
  ggtitle("Enrichment of genic consensus peaks within ChIP/eCLiP sites")





rbind(mutate(dd3, group = "Genic"), mutate(dd2, group = "Intergenic")) %>% 
  pivot_longer(cols = contains("-")) %>% 
  mutate(name2 = str_c(name, group, sep = " - ")) %>% 
  filter(ip == "eCLiP") %>% 
  ggplot(aes(x = name2, y = value, color = group)) +
  geom_jitter(width = .15) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("S9.6-only - Genic", "S9.6-only - Intergenic")
    )
  ) +
  theme_gray(14) +
  ggpubr::rotate_x_text(30) +
  ylab("Enrichment (log2 odds ratio)") +
  ggtitle("Enrichment of consensus peaks within eCLiP sites")

rbind(mutate(dd3, group = "Genic"), mutate(dd2, group = "Intergenic")) %>% 
  pivot_longer(cols = contains("-")) %>% 
  mutate(name2 = str_c(name, group, sep = " - ")) %>% 
  # filter(ip == "ChIP") %>% 
  ggplot(aes(x = name2, y = value, color = group)) +
  geom_violin(trim = F) +
  geom_jitter(width = .15) +
  facet_wrap(~ip) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("S9.6-only - Genic", "S9.6-only - Intergenic")
    )
  ) +
  theme_bw(14) +
  scale_y_continuous(limits=c(-.5, 8), expand = c(0,0)) +
  xlab(NULL) +
  ggpubr::rotate_x_text(30) +
  ylab("Enrichment (log2 odds ratio)") +
  ggtitle("Enrichment of consensus peaks within ChIP sites")





