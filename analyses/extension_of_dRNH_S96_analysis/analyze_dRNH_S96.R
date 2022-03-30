## This is an extension of the analysis from the original paper 
## Purpose is to characterize S96 and dRNH peaks and dig deeper into them

library(tidyverse)

#### MAGIC ####

ndrnh <- 56  # Number of dRNH passing QC filters
ns96 <- 191  # Number of S9.6 passing QC filters
annFull <- RLHub::annots_full_hg38()

###############

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
int <- valr::bed_intersect(s96_cons, drnh_cons)
cons2 <- list(
  "dRNH-only" = drnh_cons[! drnh_cons$name %in% int$name.y,],
  "S9.6-only" = s96_cons[! s96_cons$name %in% int$name.x,],
  "dRNH-shared" = drnh_cons[drnh_cons$name %in% int$name.y,],
  "S9.6-shared" = s96_cons[s96_cons$name %in% int$name.x,],
  "Shared" = int %>%
    mutate(start = ifelse(start.x < start.y, start.x, start.y),
           end = ifelse(end.x > end.y, end.x, end.y),
           plx = start.x + peak.x,
           ply = start.y + peak.y,
           pl = (plx + ply) / 2,
           peak = round(pl - start),
           name = paste0("Shared__", row_number())) %>%
    select(chrom, start, end, peak, name)
)


## Density plot of peak sizes
pltdat <- bind_rows(cons2, .id = "group") %>%
  mutate(
    width = end - start
  ) %>% 
  filter(group != "Shared")
mu <- pltdat %>% group_by(group) %>% summarise(mu=mean(width))
plt <- ggplot(pltdat, aes(x = width, color = group, fill = group)) +
  geom_density(alpha=.3, adjust=2) +
  scale_x_log10(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic(base_size = 18) +
  ylab("Frequency") + xlab("Peak width (bp)") +
  guides(fill=guide_legend(title=NULL), color=guide_legend(title=NULL))
plt


## Compare samples around genomic features
consInt2 <- lapply(cons2, function(x) {
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
  ggplot(aes(x = source, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  coord_flip() +
  ylab("Peak location (%)") +
  xlab(NULL) +
  scale_fill_manual(
    values = txfeat_cols,
    guide = guide_legend(title = "Feature", reverse = TRUE)
  )
plt


## Density plot of peak scores
minmax <- function(x) (x- min(x)) /(max(x)-min(x))
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


## Hockey plot of conservation across peaks
plt <- pltdat %>% 
  mutate(
    source2 = gsub(name, pattern = "__.+", replacement = "")
  ) %>% 
  ggplot(
    aes(
      x = pct_rank, y = pct_cons, color = group
    )
  ) +
  geom_line() +
  ggtitle("Conservation level of consensus peaks", 
          subtitle = "Grouped by source and overlap") +
  xlab("Proportion of total consensus peaks") +
  ylab("Percent conservation across samples") +
  theme_gray(
    base_size = 13
  ) +
  facet_wrap(~source2) +
  guides(fill=guide_legend(title=NULL), color=guide_legend(title=NULL))
plt


## TSS / TES enrichment metaplot

library(ChIPseeker)
pks <- pltdat %>% 
  group_by(group) %>% 
  {
    setNames(
      group_split(.),
      nm = group_keys(.)[[1]]
    )
  } %>% 
  as.list() %>% 
  lapply(GenomicRanges::makeGRangesFromDataFrame, keep.extra.columns = TRUE)

## TSS
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
pms <- lapply(
  pks, 
  function(x) {
    getTagMatrix(
      peak = x, TxDb = txdb, 
      upstream = 5000, downstream = 5000, 
      type = "start_site", by = "gene", 
      weightCol = "pct_cons",
    )
  }
)

plotAvgProf(pms, xlim=c(-5000, 5000), free_y = F,  facet = "row", origin_label = "TSS")
tagHeatmap(pms, xlim=c(-5000, 5000), color=NULL)

## TTS
tms <- lapply(
  pks, 
  function(x) {
    getTagMatrix(
      peak = x, TxDb = txdb, 
      upstream = 5000, downstream = 5000, 
      type = "end_site", by = "gene", 
      weightCol = "pct_cons",
    )
  }
)

plotAvgProf(tms, xlim=c(-5000, 5000), free_y = F, facet = "row", origin_label = "TTS")
tagHeatmap(tms, xlim=c(-5000, 5000), color=NULL)

## Body
bms <- lapply(
  pks, 
  function(x) {
    getTagMatrix(
      peak = x, TxDb = txdb, 
      upstream = 5000, downstream = 5000, 
      type = "end_site", by = "intron", 
      weightCol = "pct_cons"
    )
  }
)

plotAvgProf(bms, xlim=c(-5000, 5000), free_y = F, facet = "row", origin_label = "Intron")


### Genomic coverage plot
peakAnnoList <- lapply(pks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=T)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)


### Coverage
pks2 <- lapply(pks, GenomeInfoDb::keepStandardChromosomes, pruning.mode="coarse")
cplt <- covplot(pks2, weightCol="pct_cons")
cplt$layers[[1]]$aes_params = list(alpha=0)
cplt



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
db <- "GO_Biological_Process_2021"
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
      combined_score > 300 ~ 300,
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
  filter(! is.na(Term)) %>%
  ggplot(aes(x = Term, color = -log10(padj), size=combined_score, y = quant)) +
  geom_point() +
  # geom_col(position = position_dodge(.2)) +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_viridis_c(direction = -1, option = "A", end = .9,
                        guide = guide_colorbar(title = "P adj (-log10)")) 
plt









