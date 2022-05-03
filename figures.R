### Main script for generating figures (eventually should be notebook) ###
library(reticulate)
library(tidyverse)
library(kableExtra)
library(tidymodels)
library(GetoptLong)
library(enrichR)
library(VennDiagram)
library(RLSeq)
library(plotly)
library(ggprism)
source("utils.R")
# Use venv
reticulate::use_virtualenv("venv/")

# Load data from RLHub
sampleAnno <- RLHub::feat_enrich_samples()
rlrAnno <- RLHub::feat_enrich_rlregions()
rlsamples <- RLHub::rlbase_samples()
rlbps <- RLHub::rlbps()
rlr <- RLSeq::RLRangesFromRLBase("SRX1070676")
rlr <- RLSeq::featureEnrich(rlr, annotype = "full")
rlfsRes <- RLHub::rlfs_res()
rlregions <- RLHub::rlregions_meta()
annFull <- RLHub::annots_full_hg38()
cts <- RLHub::rlregions_counts()

## Genes
genes <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, 
  keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86), 
  columns = c("SEQNAME", "GENESEQSTART", "GENESEQEND", "SYMBOL")
) %>% dplyr::rename(chrom = SEQNAME, start = GENESEQSTART, end = GENESEQEND) %>%
  filter(! grepl(GENEID, pattern = "LRG.+")) %>%
  select(-GENEID) %>% mutate(chrom = paste0("chr", chrom)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# Create output directories
dir.create("results", showWarnings = FALSE)

# Human genes
humanGenes <- AnnotationDbi::select(
  x = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  columns = "SYMBOL",
  keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
)

# Block logs
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

### Figure 1/S1 ###

## Table S1 -- univariate analysis ##
rlsamples %>%
  tableone::CreateTableOne(
    data = .,
    test = FALSE, 
    vars = c("genome", "mode", "ip_type", "label")
  ) %>%
  tableone::kableone() %>% 
  kableExtra::kable_styling()  # Copy into Table S1

## Table S1 -- RLBase samples
dir.create("results/Table_S1", showWarnings = FALSE)
rlsamples %>%
  select(
    `SRA Sample ID` = rlsample,
    Prediction=prediction, 
    Label=label,
    `# of peaks` = numPeaks,
    Mode=mode,
    `IP Type`=ip_type,
    Condition=condition,
    Tissue=tissue,
    Genotype=genotype,
    Other=other,
    PMID,
    `Control Sample ID`=control,
    `Read Length` = read_length,
    `Paired End` = paired_end,
    Stranded = strand_specific,
    `Discarded during model building` = discarded,
    `Peaks URL` = peaks_s3,
    `Coverage URL` = coverage_s3,
    `RLSeq Report URL` = report_html_s3
  ) %>%
  write_csv(file = "results/Table_S1/sample_catalog.csv")

## Fig. 1C RLFS plots
dats <- tribble(
  ~rlsample, ~col,
  "SRX1025894", "#43abe1", 
  "SRX1025896", "#79919b",
  "SRX2683605", "#e662a4",
  "SRX2675009", "#a08a96"
)
rlfsTbl <- lapply(
  seq(dats$rlsample), function(i) {
    samp <- dats$rlsample[i]
    rlr <- RLRangesFromRLBase(samp)
    zscore <-  rlr@metadata$results@rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores
    tibble(
      Zscore = (zscore - min(zscore))/(max(zscore) - min(zscore)),
      shifts = rlr@metadata$results@rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifts,
      rlsample = samp
    )
  }
) %>% bind_rows()
toplt <- rlfsTbl %>%
  inner_join(rlsamples, by = "rlsample") 
plt <- toplt %>%
  dplyr::rename(Label = label, Mode = mode) %>%
  mutate(Label = factor(Label, levels = c("POS", "NEG")),
         rlsample = factor(rlsample, levels = dats$rlsample)) %>%
  ggplot(aes(x = shifts, y = Zscore, color = rlsample)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  geom_line(size = .6) +
  # facet_grid(Mode ~ Label, labeller = label_both) +
  facet_wrap(~rlsample, nrow = 1) +
  ylab("Z-score (scaled)") +
  xlab("Distance from RLFS (bp)") +
  theme_prism(border = TRUE, base_size = 10) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_manual(values = setNames(
    dats$col, nm = dats$rlsample
  )) +
  theme(legend.position = "none") +
  ggtitle("RLFS Z-Score Distribution")
plt
ggsave(plt, filename = "results/Figure_1/rlfs_zscore_plot__examples.svg", height = 4, width = 13)


## Get the summary of the datasets ##
agsmall <- RLSeq:::available_genomes %>% dplyr::select(genome=UCSC_orgID, Organism=organism)
rlsamplesNow <- rlsamples
modeDat <-  dplyr::mutate(rlsamplesNow, Mode = ifelse(mode %in% RLSeq:::auxdata$mode_cols$mode, mode, "misc")) %>%
  group_by(Mode) %>% tally()
genDat <- dplyr::rename(rlsamplesNow, Genome = genome) %>% group_by(Genome) %>% tally()
labelDat <- dplyr::rename(rlsamplesNow, Label = label) %>% group_by(Label) %>% tally()
predDat <- dplyr::rename(rlsamplesNow, Prediction = prediction) %>% 
  group_by(Prediction) %>%
  tally() %>%
  mutate(Prediction = ifelse(is.na(Prediction), "No Peaks", Prediction))
ipDat <- dplyr::rename(rlsamplesNow, IP_type = ip_type) %>% 
  mutate(IP_type = ifelse(IP_type %in% c("None", "RNA-m6A"), "Other", IP_type)) %>%
  group_by(IP_type) %>%
  tally()
datList <- list("Mode"=modeDat, "Label"=labelDat, "Prediction"=predDat, "Genome"=genDat, "IP_type" = ipDat)
outpath <- "results/Figure_1/"
auxdata <- RLSeq:::auxdata
ip_cols <- c("dRNH" = "#c386c4", "S9.6"="#82d0e8")

auxdata$ip_type_cols <- tribble(
  ~ip_type, ~col,
  "dRNH", "#c386c4",
  "S9.6", "#82d0e8",
  "Other", "#e0dede"
)
dir.create(outpath, showWarnings = FALSE)
pltLst <- lapply(names(datList), function(dat) {
  dataNow <- datList[[dat]]
  if (dat != "Genome") {
    cs <- auxdata[[paste0(tolower(dat), "_cols")]]
    mark <- list(colors=setNames(cs[,'col', drop=TRUE], nm=cs[,tolower(dat),drop=TRUE])[dataNow[,dat,drop=T]])
  } else {
    mark <- NULL
  }
  plt <- plot_ly(type = "pie") %>%
    add_pie(data=dataNow, labels = dataNow[,dat,drop=T], values = ~n, textinfo='label+value',
            marker = mark, 
            insidetextorientation='horizontal', hole=.6, rotation=250) %>%
    layout(showlegend = FALSE, title=list(text = dat, x=0.15), margin = list(l = 100, r = 100, t=100, b=100),
           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  save_image(plt, file = file.path(outpath, paste0("donut__", dat, ".svg")), format = "svg")
  
})
names(pltLst) <- names(datList)


# Num peaks across conditions -- pertains to Figure S1
pltdat <- rlsamples %>%
  filter(label == "POS") %>%
  mutate(numPeaks = numPeaks / 1000) %>%
  group_by(mode) %>%
  mutate(numpeakmed = median(numPeaks, na.rm = TRUE)) %>%
  ungroup()
ord <- pltdat %>%
  group_by(mode) %>%
  summarise(npeak = mean(numPeaks))
plt <- pltdat %>%
  filter(! is.na(numPeaks), numPeaks > 0) %>%
  arrange(numpeakmed) %>%
  mutate(mode = factor(mode, levels = unique(mode))) %>%
  ggplot(mapping = aes(x = mode, fill = mode,
                       y = numPeaks)) +
  geom_boxplot(color = "black", width = .7, outlier.shape=NA) +
  geom_jitter(size=1, width=.2, alpha=.75) +
  ylab("Peaks Called (thousands)") +
  xlab(NULL) +
  labs(title = "Peaks called by mode") +
  scale_fill_manual(values = setNames(RLSeq:::auxdata$mode_cols$col, nm = RLSeq:::auxdata$mode_cols$mode)) +
  guides(y = "prism_offset_minor") + 
  theme_prism(base_size = 14) + 
  ggpubr::rremove("legend") +
  theme(axis.text.x = element_text(angle = 36, vjust = 1,
                                   hjust = 1)) 



ggsave(plt, filename = "results/Figure_1/peaks_called_by_mode.pdf", height = 3.5, width = 12)


## Figure showing sizes by modality
fls <- list.files(
  "../RLBase-data/rlbase-data/rlpipes-out/peaks/", 
  pattern = ".broadPeak$", full.names = TRUE
)
pkwd <- parallel::mclapply(
  seq(nrow(pltdat)), function(i) {
    message(i)
    samp <- pltdat$rlsample[i]
    pk <- fls[grep(fls, pattern = samp)]
    if (length(pk) == 0) return(NULL)
    e <- try(
      {
        valr::read_broadpeak(pk[1]) %>% 
          mutate(width = end - start) %>% 
          summarise(
            medwidth = median(width),
            meanwidth = mean(width),
            medqval = median(qvalue),
            meanqval = mean(qvalue),
            corwidqval_p = cor.test(width, qvalue, method = "spearman")$p.value,
            corwidqval_r = cor.test(width, qvalue, method = "spearman")$estimate
          ) %>% 
          mutate(rlsample = samp)
      }
    )
    if (class(e)[1] != "try-error") return(e) else return(NULL)
  }, mc.cores = 44
) %>% bind_rows()

samp <- pkwd %>%filter(corwidqval_p == 0) %>% 
  slice_max(corwidqval_r) 
plt1 <- valr::read_broadpeak(fls[grep(fls, pattern = samp$rlsample)]) %>% 
  mutate(width = end - start) %>%
  ggplot(aes(x = width, y = qvalue)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10() +
  theme_gray(base_size = 10) +
  xlab("Peak width (log10)") +
  ylab("Peak Q-Value (log10)") +
  ggtitle(
    paste0("Peak width vs strength (", samp$rlsample, ")"),
    subtitle = paste0("Rho: ", round(samp$corwidqval_r, 4), " -- padj: ", samp$corwidqval_p)
  )
plt2 <- pkwd %>% 
  mutate(
    padj = p.adjust(corwidqval_p),
    sig = padj < 0.05
  ) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(pct = 100*(n/sum(n))) %>% 
  ggplot(aes(x=1, y = pct, fill=sig)) +
  geom_bar(stat = "identity", width=.4, position = "stack") +
  coord_flip() +
  xlab(NULL) +
  ylab("Percentage of samples") +
  theme_gray(base_size = 10) +
  scale_fill_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "grey")
  ) +
  guides(fill = guide_legend(title = "Significant\ncorrelation")) +
  ggtitle("Pct of samples w/ significant correlation \nbetween width and strength of peaks")

# Make plot
pltdat2 <- inner_join(pkwd, pltdat, by = "rlsample") %>%
  group_by(mode) %>%
  mutate(
    pkwdmed = median(medwidth, na.rm=TRUE)
  )
plt3 <- pltdat2 %>%
  filter(! is.na(medwidth), medwidth > 0, ! is.na(pkwdmed), label == "POS") %>%
  arrange(pkwdmed) %>%
  mutate(mode = factor(mode, levels = unique(.$mode))) %>%
  ggplot(mapping = aes(x = mode, fill = mode,
                       y = medwidth)) +
  geom_boxplot(color = "black", width = .7, outlier.shape=NA) +
  geom_jitter(size=1, width=.2, alpha=.75) +
  ylab("Median peak width (base pairs)") +
  xlab(NULL) +
  labs(title = "Median peak widths by mode") +
  scale_fill_manual(values = setNames(RLSeq:::auxdata$mode_cols$col, nm = RLSeq:::auxdata$mode_cols$mode)) +
  guides(y = "prism_offset_minor") + 
  theme_prism(base_size = 14) + 
  ggpubr::rremove("legend") +
  theme(axis.text.x = element_text(angle = 36, vjust = 1,
                                   hjust = 1)) 

ggpubr::ggarrange(
  plt1, plt2, plt3, labels = "AUTO", widths = c(2, 2, 3), nrow = 1
)

### Figure 2 / S2 / S3 ###
dir.create("results/Figure_2/", showWarnings = FALSE)

# predlab colors
predlabcond <- tribble(
  ~cond, ~prediction,  ~label,
  "POS_NEG", "POS", "NEG",
  "POS_POS", "POS", "POS",
  "NEG_POS", "NEG", "POS",
  "NEG_NEG", "NEG", "NEG"
)
predlabcols <- tribble(
  ~cond, ~col,
  "NEG_NEG", "#8a2c2c",
  "NEG_POS", "#A76767",
  "POS_NEG", "#7aa4c4",
  "POS_POS", "#2270ab"
)

## Metaplots of RLFS analysis (Supplemental)
rlfsTbl <- lapply(
  seq(rlfsRes), function(i) {
    rlfsNow <- rlfsRes[[i]]
    if (is.na(rlfsNow)) return(NULL)
    zscore <- rlfsNow$rlfsData$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores
    # Min max scaling to 0-1 range
    tibble(
      Zscore = (zscore - min(zscore))/(max(zscore) - min(zscore)),
      shifts = rlfsNow$rlfsData$`Z-scores`$`regioneR::numOverlaps`$shifts,
      rlsample = names(rlfsRes)[i]
    )
  }
)
rlfsTbl <- rlfsTbl[-sapply(rlfsTbl, is.null)] %>% bind_rows()
toplt <- rlfsTbl %>%
  inner_join(rlsamples, by = "rlsample") 

## Fig S3 ##

# With just selected samples
pos_pos <- c(
  "SRX3084731",
  "SRX4776642",
  "SRX1070677",
  "SRX8961684",
  "SRX2683608",
  "SRX7347914"
)
pos_neg <- c(
  "SRX10229652",
  "SRX5696405",
  "SRX3892923",
  "ERX3974986",
  "ERX3974965",
  "SRX6427715"
)
neg_pos <- c(
  "SRX534096",
  "SRX6686232",
  "SRX8110930",
  "SRX4732954",
  "SRX2733012",
  "SRX2642935"
)
neg_neg <- c(
  "SRX2187015",
  "SRX7805788",
  "SRX1029478",
  "SRX1070681",
  "SRX6779956",
  "SRX2187019"
)
samps <- c(pos_pos, pos_neg, neg_pos, neg_neg)
plt1 <- toplt %>%
  filter(
    rlsample %in% samps
  ) %>%
  inner_join(predlabcond) %>%
  arrange(cond) %>%
  mutate(rlsample = factor(rlsample, levels = unique(rlsample))) %>%
  ggplot(aes(x = shifts, y = Zscore, color = cond, fill=cond)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  geom_line() +
  facet_wrap(rlsample~., nrow = 4, ncol = 6) +
  ylab("Z-score (scaled)") +
  xlab("Distance from RLFS (bp)") +
  theme_prism(border = TRUE, base_size = 10) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  scale_fill_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  theme(axis.text.x = element_text(size = 8)) + ggpubr::rremove("legend")
plt2 <-  toplt %>%
  filter(
    rlsample %in% samps
  ) %>%
  inner_join(predlabcond) %>%
  arrange(cond) %>%
  mutate(rlsample = factor(rlsample, levels = unique(rlsample))) %>%
  ggplot(aes(x = shifts, y = Zscore, color = cond, fill=cond, group=rlsample)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  geom_line() +
  facet_wrap(cond~., nrow = 4, ncol = 1) +
  ylab(NULL) +
  xlab(" ") +
  theme_prism(border = TRUE, base_size = 10) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  scale_fill_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  ggpubr::rremove("legend") +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

# Meta
plt3 <- toplt %>%
  inner_join(predlabcond) %>%
  ggplot(aes(x = shifts, y = Zscore, color = cond, fill=cond)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  stat_smooth(method="loess", span=.01, se=TRUE, alpha=.3) +
  facet_wrap(cond~., nrow = 4, ncol = 1) +
  ylab(NULL) +
  xlab(" ") +
  theme_prism(border = TRUE, base_size = 10) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  scale_fill_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  theme(axis.text.x = element_text(size = 8), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ggpubr::rremove("legend")

plt <- ggpubr::ggarrange(plt1, plt2, plt3, ncol = 3, widths = c(6, 1, 1))
ggsave(plt, filename = "results/Figure_2/rlfs_zscore_supplement.svg", height = 8, width = 16)


## Figure 2 ##

# Metaplots of RLFS (main panel)
plt <- toplt %>%
  inner_join(predlabcond) %>%
  rename(Label = label, Prediction=prediction) %>%
  ggplot(aes(x = shifts, y = Zscore, color = cond, fill=cond)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  stat_smooth(method="loess", span=.01, se=TRUE, alpha=.3) +
  facet_grid(Prediction ~ Label, labeller = label_both) +
  ylab("Z-score (scaled)") +
  xlab("Distance from RLFS (bp)") +
  theme_prism(border = TRUE, base_size = 10) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  scale_fill_manual(values = setNames(predlabcols$col, nm = predlabcols$cond)) +
  theme(legend.position = "none") +
  ggtitle("RLFS Z-Score Distribution", "Label vs Prediction")
ggsave(plt, filename = "results/Figure_2/rlfs_zscore_plot.svg", height = 5, width = 7.5)

rlsamples %>%
  filter(! is.na(prediction)) %>%
  group_by(label, prediction) %>% 
  tally() %>%
  ungroup() %>%
  mutate(pct = n/sum(n))

## Figure S2 ##

# Percentage label and prediction
colors1<- colorRampPalette( RColorBrewer::brewer.pal(9, "Blues"))(255)[1:150]
annorow <- rlsamples %>%
  group_by(mode) %>% 
  summarise(`Total (n)`=n()) %>%
  column_to_rownames(var = "mode")
annocol <- tribble(
  ~row, ~label, ~prediction,
  "POS_NEG", "POS", "NEG",
  "POS_POS", "POS", "POS",
  "NEG_POS", "NEG", "POS",
  "NEG_NEG", "NEG", "NEG"
) %>%
  column_to_rownames("row")
heatdata <- rlsamples %>%
  group_by(mode, label, prediction) %>%
  tally() %>%
  na.omit() %>%
  ungroup() %>%
  group_by(mode) %>%
  mutate(pct = (n / sum(n)) * 100) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(label, prediction), names_from = mode, values_from = pct, values_fill = 0) %>%
  mutate(row = paste0(label, "_", prediction)) %>%
  arrange(desc(label), desc(prediction)) %>%
  select(-label, -prediction) %>%
  column_to_rownames(var = "row") %>%
  t() 
cat_cols <- list(
  "label" = setNames(RLSeq:::auxdata$label_cols$col, nm = RLSeq:::auxdata$label_cols$label),
  "prediction" = setNames(RLSeq:::auxdata$prediction_cols$col, nm = RLSeq:::auxdata$prediction_cols$prediction),
  "Total (n)" = RColorBrewer::brewer.pal(n = 9, "OrRd")
)
chm <- ComplexHeatmap::pheatmap(
  heatdata, color = colors1, main = "Percentage samples labels by Mode",
  fontsize = 13, cluster_rows = FALSE, cluster_cols = FALSE,
  display_numbers = TRUE, fontsize_number = 10.5,
  annotation_col = annocol[colnames(heatdata), , drop=FALSE], 
  show_colnames = FALSE, name = "Representation (%)",
  annotation_colors = cat_cols, use_raster=FALSE,
  annotation_row = annorow[rownames(heatdata),,drop=FALSE],
  number_color = "black", angle_col = "0"
)

pdf(qq("results/Figure_2/heatmap_representation.pdf"), 
    width = 8, height = 8)
ComplexHeatmap::draw(chm)
dev.off()



### Figure 3/S4
dir.create("results/Figure_3", showWarnings = FALSE)

## Genome Browser -- get scaling factors for bigWigs
## From https://www.biostars.org/p/413626/#414440
rlcts <- RLHub::rlregions_counts()
cts_mat <- rlcts@assays@data$cts
# normFact <- edgeR::calcNormFactors(object = cts_mat, method="TMM")
# libSize <- colSums(cts_mat)
# sizeFactors <- normFact * libSize / 1000000
# sizeFactors.Reciprocal <- 1/sizeFactors
# sizeFactTbl <- tibble(
#   "rlsample" = names(sizeFactors),
#   "sizeFactors" = sizeFactors,
#   "sizeFactors.Reciprocal" = sizeFactors.Reciprocal
# ) %>% write_csv(file = "data/sizeFactors.csv")

########## DOES NOT NEED TO BE RERUN WHEN NEW DATA ADDED ##############

# Re-run deepTools with the scaling (done from the CLI)
# 1. Predict fragment length
# macs3 predictd -i "../RLBase-data/rlbase-data/rlpipes-out/bam/SRX2455193/SRX2455193_hg38.bam" -g hs
# tribble(
#   ~rlsample, ~predict_fragsize,
#   "ERX3974964", 245,
#   "ERX3974965", 255,
#   "SRX2455189", 197,
#   "SRX2455193", 182,
#   "SRX1025894", 191,
#   "SRX1025896", 294,
#   "SRX2683605", 296,
#   "SRX2675009", 260
# )
# 2. Run deeptools with scaling factors and extended fragments
# bamCoverage -b "bam/ERX3974964/ERX3974964_hg38.bam" -o coverage_scaled/ERX3974964_hg38.scale.bw --scaleFactor 0.07312689 -p 44 -e 245 --ignoreDuplicates
# bamCoverage -b "bam/ERX3974965/ERX3974965_hg38.bam" -o coverage_scaled/ERX3974965 --scaleFactor 0.2667499 -p 44 -e 255 --ignoreDuplicates
# bamCoverage -b "bam/SRX2455189/SRX2455189_hg38.bam" -o coverage_scaled/SRX2455189_hg38.scaled.bw --scaleFactor 0.7608039 -p 44 -e 197 --ignoreDuplicates
# bamCoverage -b "bam/SRX2455193/SRX2455193_hg38.bam" -o coverage_scaled/SRX2455193_hg38.scaled.bw --scaleFactor 0.5406119 -p 44 -e 182 --ignoreDuplicates
# bamCoverage -b "bam/SRX1025894/SRX1025894_hg38.bam" -o coverage_scaled/SRX1025894_hg38.scale.bw --scaleFactor 0.1704441 -p 44 -e 191 --ignoreDuplicates
# bamCoverage -b "bam/SRX1025896/SRX1025896_hg38.bam" -o coverage_scaled/SRX1025896_hg38.scale.bw --scaleFactor 0.7183404 -p 44 -e 294 --ignoreDuplicates
# bamCoverage -b "bam/SRX2683605/SRX2683605_hg38.bam" -o coverage_scaled/SRX2683605_hg38.scale.bw --scaleFactor 0.6507813 -p 44 -e 296 --ignoreDuplicates
# bamCoverage -b "bam/SRX2675009/SRX2675009_hg38.bam" -o coverage_scaled/SRX2675009_hg38.scale.bw --scaleFactor 0.9435923 -p 44 -e 260 --ignoreDuplicates


##########################################################


# Gold standard heatmap
heatdata <- RLSeq::corrHeatmap(rlr, returnData = TRUE)
hm <- ComplexHeatmap::Heatmap(
  heatdata$corrRes,
  col = heatdata$continuous_pal,
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  use_raster = TRUE,
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    df = heatdata$annoCorr %>% select(-Selected), col = heatdata$cat_cols
  ),
  name = "Corr (R)"
)
pdf(qq("results/Figure_3/corr_heatmap.pdf"), 
    width = 11, height = 8)
ComplexHeatmap::draw(hm)
dev.off()

## Show the strength of these results with the top hit
pltdat <- RLSeq::plotEnrichment(
  rlr, returnData = TRUE, 
  pred_POS_only = FALSE,
  rlbaseRes = sampleAnno, 
  rlsamples = rlsamples,
  label_POS_only = FALSE
) 

# Genes
plt <- plot_multi_feature(
  features = unique(pltdat$Transcript_Features$type),
  db = "Transcript_Features", 
  factororder = c("TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS"),
  pltdat = pltdat,
  lmts = c(-6, 12),
  axislmts = c(-6, 15),
  yb = 11.5
) +
  ggtitle("Genic Feature Enrichment in RLBase samples", subtitle = NULL) + xlab(NULL)
ggsave(plt, filename = "results/Figure_3/genic_enrichment.svg", height = 5, width = 11)

# G4Q
pltdat2 <- pltdat
pltdat2$G4Qpred$type <- "G4 Predicted"
pltdat2$G4Qexp$type <- "G4 ChIP"
feats <- c("G4 Predicted", "G4 ChIP")
plt <- plot_multi_feature(
  features = feats, db = c("G4Qpred", "G4Qexp"), 
  factororder = feats,
  pltdat = pltdat2,
  lmts = c(-6, 7),
  axislmts = c(-6, 8),
  yb = 11.5
) +
  ggtitle("G4Q Enrichment in RLBase samples", subtitle = NULL) + xlab(NULL)
plt
ggsave(plt, filename = "results/Figure_3/g4q_enrichment.svg", height = 5, width = 7)


# GC Skew
plt <- plot_multi_feature(
  features = c("G SKEW", "C SKEW", "CpG Islands"), db = c("skewr", "CpG_Islands"), 
  factororder = c("CpG Islands", "G SKEW", "C SKEW"),
  pltdat = pltdat,
  lmts = c(-6, 8),
  axislmts = c(-6, 10),
  yb = 11.5
) +
  ggtitle("G/C Skew Enrichment in RLBase samples", subtitle = NULL) + xlab(NULL)
plt
ggsave(plt, filename = "results/Figure_3/gc_enrichment.svg", height = 5, width = 9)

to_collect <- c(
  "SRX1070678", 
  "SRX9182646",
  "SRX8122769",
  "SRX2683605",
  "SRX5547605",
  "SRX10229647"
)
rlsamples %>% 
  filter(rlsample %in% to_collect) %>% 
  pull(peaks_s3) %>% 
  paste0(RLSeq:::RLBASE_URL, "/", .) %>% 
  cat(sep = "\n")

# track type=broadPeak name=" " useScore=0
# https://rlbase-data.s3.amazonaws.com/peaks/SRX1070678_hg38.broadPeak
# 
# track type=broadPeak name="  " useScore=0
# https://rlbase-data.s3.amazonaws.com/peaks/SRX8122769_hg38.broadPeak
# 
# track type=broadPeak name="   " useScore=0
# https://rlbase-data.s3.amazonaws.com/peaks/SRX9182646_hg38.broadPeak
# 
# track type=broadPeak name="    " useScore=0
# https://rlbase-data.s3.amazonaws.com/peaks/SRX10229647_hg38.broadPeak
# 
# track type=broadPeak name="     " useScore=0
# https://rlbase-data.s3.amazonaws.com/peaks/SRX2683605_hg38.broadPeak
# 
# track type=broadPeak name="      " useScore=0
# https://rlbase-data.s3.amazonaws.com/peaks/SRX5547605_hg38.broadPeak


### Figure 4/S5 ###
### RL Regions
dir.create("results/Figure_4", showWarnings = FALSE)

dataNow <- rlsamples %>%
  filter(
    genome == "hg38",
    label == "POS",
    ! is.na(prediction),
    ip_type %in% c("S9.6", "dRNH")
  ) %>%
  group_by(ip_type) %>%
  {setNames(group_split(.), group_keys(.)[[1]])} %>%  
  lapply(
    function(x) {
      x %>% 
        mutate(Included = ifelse(
          label == "POS" & prediction == "POS" & numPeaks > 5000, "Yes", "No"
        )) %>%
        group_by(Included) %>%
        tally()
    }
  )

ip_cols <- list("dRNH"=rev(c("Yes" = "#e2a3e3", "No" = "#d6ced6")),
                "S9.6" = rev(c("Yes"="#82d0e8", "No" = "#c7d1d4")))

plt1 <- plot_ly(type = "pie") %>%
  add_pie(data=dataNow$dRNH, labels = dataNow$dRNH$Included, values = ~n, textinfo='label+value',
          marker = list(colors=ip_cols$dRNH), 
          insidetextorientation='horizontal', hole=.6) %>%
  layout(showlegend = FALSE, title=list(text = "Included", x=0.15), margin = list(l = 100, r = 100, t=100, b=100),
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
save_image(plt1, file = "results/Figure_4/dRNH__sample_donut.svg", format = "svg")
plt2 <- plot_ly(type = "pie") %>%
  add_pie(data=dataNow$S9.6, labels = dataNow$S9.6$Included, values = ~n, textinfo='label+value',
          marker = list(colors=ip_cols$S9.6), 
          insidetextorientation='horizontal', hole=.6) %>%
  layout(showlegend = FALSE, title=list(text = "Included", x=0.15), margin = list(l = 100, r = 100, t=100, b=100),
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
save_image(plt2, file = "results/Figure_4/S96__sample_donut.svg", format = "svg")


## R-loop venn diagram
drnh_cons <- read_tsv(file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_dRNH.narrowPeak"),
                      col_names = c("chrom", "start", "end", "name", "score", "strand", "signalVal", "pVal", "qVal", "peak"), 
                      show_col_types = FALSE, progress = FALSE)
s96_cons <- read_tsv(file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_S96.narrowPeak"),
                     col_names = c("chrom", "start", "end", "name", "score", "strand", "signalVal", "pVal", "qVal", "peak"), 
                     show_col_types = FALSE, progress = FALSE)
genome <- read_tsv("https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", 
                   col_names = c("chrom", "size"), show_col_types = FALSE, progress = FALSE)
print(valr::bed_fisher(drnh_cons, y = s96_cons, genome = genome))
intcons <- valr::bed_intersect(x = drnh_cons, y = s96_cons) 
# Get size of overlap for venn diagram
c(
  GenomicRanges::makeGRangesFromDataFrame(
    data.frame(chrom=intcons$chrom, start=intcons$start.x, end=intcons$end.x)
  ),
  GenomicRanges::makeGRangesFromDataFrame(
    data.frame(chrom=intcons$chrom, start=intcons$start.y, end=intcons$end.y)
  )
) %>% GenomicRanges::reduce() %>% length()  # Shared number
# Get size of difference
length(which(! drnh_cons$name %in% intcons$name.x))  # dRNH-unique number
length(which(! s96_cons$name %in% intcons$name.y))  # S96-unique number


## Density plot of peak sizes
pltdat <- bind_rows(list(
  tibble(ip_type = "dRNH", width = drnh_cons$end-drnh_cons$start),
  tibble(ip_type = "S9.6", width = s96_cons$end-s96_cons$start)
)) %>%
  filter(width < 80000) 
mu <- plyr::ddply(pltdat, "ip_type", summarise, grp.mean=mean(width))
plt <- ggplot(pltdat, aes(x = width, color = ip_type, fill = ip_type)) +
  geom_density(alpha=.3, adjust=2) +
  scale_x_log10(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic(base_size = 18) +
  scale_fill_manual(values = c("dRNH" = "#e2a3e3", "S9.6" = "#82d0e8")) +
  scale_color_manual(values = c("dRNH" = "#e2a3e3", "S9.6" = "#82d0e8")) +
  ylab("Frequency") + xlab("Peak width (bp)") +
  guides(fill=guide_legend(title=NULL), color=guide_legend(title=NULL))
ggsave(plt, filename = "results/Figure_4/width_plot.svg")

## Metaplot of coverage around genes for S9.6 and dRNH coverage
# Need to run the line below to generate signal matrix
# computeMatrix scale-regions -S ../RLBase-data/rlbase-data/rlregions/rlregions_S96.bw ../RLBase-data/rlbase-data/rlregions/rlregions_dRNH.bw -R ~/.rlpipes_genomes/hg38/hg38.ensGene.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 -o ../RLBase-data/misc-data/dRNH_S96_matrix.mat.gz -p 44 --verbose
metap <- read_tsv("../RLBase-data/misc-data/dRNH_S96_matrix.mat.gz", skip = 1, col_names = FALSE)
metadat <- read_lines("../RLBase-data/misc-data/dRNH_S96_matrix.mat.gz", n_max = 1)
metadat <- gsub(metadat, pattern = "@", replacement = "") %>% jsonlite::fromJSON()
ndrnh <- dataNow$dRNH$n[dataNow$dRNH$Included == "Yes"]
ns96 <- dataNow$S9.6$n[dataNow$S9.6$Included == "Yes"]
mat <- metap %>%
  select(-c(1:6)) 
mat %>%
  colMeans(na.rm = TRUE) -> cm
divs <- setNames(c(
  ndrnh,
  ns96
), nm = c("rlregions_dRNH", "rlregions_S96"))
pltdats <- tibble("score" = cm,
                  "loc" = rep(seq(1100), 2),
                  "div" = c(rep(divs[metadat$sample_labels[1]], metadat$sample_boundaries[2]),
                            rep(divs[metadat$sample_labels[2]], metadat$sample_boundaries[2])),
                  "label" = c(rep(metadat$sample_labels[1], metadat$sample_boundaries[2]),
                              rep(metadat$sample_labels[2], metadat$sample_boundaries[2]))) %>%
  mutate(
    score = score / div,
    score = score / sum(score),
    label = factor(label, levels = rev(unique(label)))) %>%
  arrange(desc(label)) 
mcm <- min(pltdats$score)
mxcm <- max(pltdats$score)
mcdiff <- mxcm-mcm
plt <- pltdats %>%
  ggplot(aes(x = loc, y = score, color = label, fill=label)) +
  geom_vline(xintercept = 300,  color="grey") +
  geom_vline(xintercept = 800,  color="grey") +
  geom_line() +
  theme_classic(base_size = 14) +
  ylab("Density") +
  xlab("Gene position (bp)") +
  scale_x_continuous(
    breaks = c(1, 300, 800, 1100),
    labels = c("-3000", "TSS", "TTS", "+3000"),
    expand = c(0,0)
  ) +
  ggpubr::rremove("legend") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(mcm, mxcm),
                     breaks = mcm + c(0*mcdiff, .25*mcdiff, .5*mcdiff, .75*mcdiff, 1*mcdiff),
                     labels = c(0, .25, .5, .75, 1)) 
ggsave(plt, filename = "results/Figure_4/gene_position_rlregions.svg", height = 4.35, width = 6.0)


## Compare samples around genomic features
cons <- list(
  "dRNH" = drnh_cons,
  "S9.6" = s96_cons
)
consInt <- lapply(cons, function(x) {
  if (x$name[1] == "dRNH__1") {
    x <- mutate(
      x, 
      pct_cons = 100 * ((score / 10) / ndrnh),  
      source = "dRNH"
    )
  } else {
    x <- mutate(
      x, 
      pct_cons = 100 * ((score / 10) / ns96),  # Compared to 184 for S9.6
      source = "S9.6"
    )
  }
  
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
}) %>% bind_rows()

# Same but with shared and not shared
int <- valr::bed_intersect(s96_cons, drnh_cons)
cons2 <- list(
  "dRNH-only" = drnh_cons[! drnh_cons$name %in% int$name.y,],
  "S9.6-only" = s96_cons[! s96_cons$name %in% int$name.x,],
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
}) %>% bind_rows()

# Annotate TTS, TSS, and Gene Body
txfeat_cols <- setNames(
  rev(c("#b36262", "#c4cc6a", "#d17bdb", "#4bc6d6", "#83d647", "#9494e3", "#7d7d7d")),
  nm = rev(c("TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "Intergenic"))
)
txfeats <- annFull[sapply(names(annFull), function(x) grepl(x, pattern = "Transcript_Features.*"))]
txfeats <- lapply(names(txfeats), function(txname) {
  x <- txfeats[[txname]]
  x$db <- gsub(txname, pattern = "(.+)__(.+)", replacement = "\\1")
  x$type <- gsub(txname, pattern = "(.+)__(.+)", replacement = "\\2")
  return(x)
}) %>% bind_rows()
oltx <- consInt %>%
  select(chrom, start, end, name, source) %>%
  valr::bed_intersect(txfeats)
oltxsum <- oltx %>%
  group_by(name.x) %>%
  summarise(
    type = ifelse("TSS" %in% type.y, "TSS", 
                  ifelse("TTS" %in% type.y, "TTS",
                         ifelse("fiveUTR" %in% type.y, "fiveUTR", 
                                ifelse("threeUTR" %in% type.y, "threeUTR", 
                                       ifelse("Exon" %in% type.y, "Exon", "Intron")))))
  ) %>%
  mutate(source = gsub(name.x, pattern = "(.+)__.+", replacement = "\\1"))

oltxres <- oltxsum %>%
  dplyr::rename(name = name.x) %>%
  full_join(consInt, by = "name") %>%
  select(name, source = source.y, type) %>%
  mutate(type = ifelse(is.na(type), "Intergenic", type)) %>%
  distinct(name, .keep_all = TRUE)
plt <- oltxres %>%
  group_by(source, type) %>%
  tally() %>%
  mutate(n = 100*n / sum(n),
         type = factor(type, levels = rev(c(
           "TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "Intergenic"
         ))))  %>%
  ggplot(aes(x = source, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  coord_flip() +
  ylab("Peak location (%)") +
  xlab(NULL) +
  scale_fill_manual(values = txfeat_cols, guide = guide_legend(title = "Feature", reverse = TRUE))
plt
ggsave(plt, filename = "results/Figure_4/consensus__txplot.svg")

# With shared
oltx2 <- consInt2 %>%
  mutate(source = gsub(name, pattern = "(.+)__.+", replacement = "\\1")) %>%
  select(chrom, start, end, name, source) %>%
  valr::bed_intersect(txfeats)
oltxsum2 <- oltx2 %>%
  group_by(name.x) %>%
  summarise(
    type = ifelse("TSS" %in% type.y, "TSS", 
                  ifelse("TTS" %in% type.y, "TTS",
                         ifelse("fiveUTR" %in% type.y, "fiveUTR", 
                                ifelse("threeUTR" %in% type.y, "threeUTR", 
                                       ifelse("Exon" %in% type.y, "Exon", "Intron")))))
  ) %>%
  mutate(source = gsub(name.x, pattern = "(.+)__.+", replacement = "\\1"))

oltxres2 <- oltxsum2 %>%
  dplyr::rename(name = name.x) %>%
  full_join(consInt2, by = "name") %>%
  mutate(source = gsub(name, pattern = "(.+)__.+", replacement = "\\1")) %>%
  select(name, source, type) %>%
  mutate(type = ifelse(is.na(type), "Intergenic", type)) %>%
  distinct(name, .keep_all = TRUE)
plt <- oltxres2 %>%
  group_by(source, type) %>%
  tally() %>%
  mutate(
    n = 100*n / sum(n),
    source = ifelse(
      source == "dRNH", "dRNH-only",
      ifelse(
        source == "S96", "S9.6-only", source
      )
    ),
    source = factor(source, levels = rev(c(
      "dRNH-only", "S9.6-only", "Shared"
    ))),
    type = factor(type, levels = rev(c(
      "TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "Intergenic"
    )))
  )  %>%
  ggplot(aes(x = source, y = n, fill = type)) +
  geom_bar(stat = "identity", color="black", position = "stack") +
  theme_bw(base_size = 20) +
  coord_flip() +
  ylab("Peak location (%)") +
  xlab(NULL) +
  ggtitle("Consensus site overlap in genes") +
  scale_fill_manual(values = txfeat_cols, guide = guide_legend(title = "Feature", reverse = TRUE))
plt
ggsave(plt, filename = "results/Figure_4/consensus__txplot_2.svg", height = 6, width = 8)

## Venn showing olap of S96, dRNH, and RLFS
# ol <- ChIPpeakAnno::findOverlapsOfPeaks(
#   GenomicRanges::makeGRangesFromDataFrame(cons$dRNH),
#   GenomicRanges::makeGRangesFromDataFrame(cons$S9.6),
#   RLSeq:::getRLFSAnno(object = rlr)
# )
# pdf("results/Figure_4/RLFS_venn.pdf", height = 4, width = 7)
# ChIPpeakAnno::makeVennDiagram(ol, fill = c("#e2a3e3", "#82d0e8", "white"), 
#                               NameOfPeaks = c("dRNH", "S9.6", "RLFS"))
# dev.off()

## Venn showing olap of dRNH and S96
ol <- ChIPpeakAnno::findOverlapsOfPeaks(
  GenomicRanges::makeGRangesFromDataFrame(cons$dRNH),
  GenomicRanges::makeGRangesFromDataFrame(cons$S9.6)
)
pdf("results/Figure_4/d96_dRNH_venn.pdf", height = 4, width = 7)
ChIPpeakAnno::makeVennDiagram(ol, fill = c("#e2a3e3", "#82d0e8"), 
                              NameOfPeaks = c("dRNH", "S9.6"))
dev.off()


### Figure 5/S6 ###

## PCA from RLRegion counts (Junk included!)
ip_correct <- rlsamples$rlsample[rlsamples$ip_type %in% c("S9.6", "dRNH")]
ctsPos <- cts[
  gsub(rownames(cts), pattern = "All_", replacement = "") %in% rlregions$rlregion[rlregions$source == "dRNH S96"]
  ,
  # ! is.na(cts$prediction) &
  # cts$prediction == "POS" & 
  cts$experiment %in% ip_correct &
    cts$label == "POS" 
  # cts$numPeaks > 5000
]
ctsmat <- ctsPos@assays@data$cts
vstmat <- DESeq2::vst(ctsmat)
pcdata <- prcomp(t(vstmat))
pcsmall <- pcdata$x[,1:7]
pctpc <- round(100 * pcdata$sdev / sum(pcdata$sdev), digits = 2)
pltdat <- pcsmall %>% 
  as.data.frame() %>%
  rownames_to_column(var = "name") %>%
  inner_join(rlsamples, by = c("name" = "rlsample")) %>%
  mutate(mode = case_when(
    mode %in% RLSeq:::auxdata$mode_cols$mode ~ mode,
    TRUE ~ "misc"
  ))
ip_cols <- c("dRNH" = "#C86EAB", "S9.6"="#0AB9E4")
pltpca <- as_tibble(pltdat) %>%
  arrange(desc(prediction)) %>% 
  filter(! is.na(prediction)) %>% 
  ggplot(aes(x = PC1, y = PC2, color = prediction, group = tissue, label = name)) +
  geom_hline(yintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
  geom_point(size=2.5) +
  facet_wrap(~ip_type, nrow = 2) +
  theme_bw(base_size = 17) +
  ggtitle("PCA with false positive samples included") +
  xlab(paste0("PC1 (", pctpc[1], "%)")) +
  ylab(paste0("PC2 (", pctpc[2], "%)")) 
pltpca
ggsave(plot = pltpca, filename = "results/Figure_5/pca_plot_with_false_pos.png", height = 10, width = 6)


## PCA from RLRegion counts
ip_correct <- rlsamples$rlsample[rlsamples$ip_type %in% c("S9.6", "dRNH")]
ctsPos <- cts[
  gsub(rownames(cts), pattern = "All_", replacement = "") %in% rlregions$rlregion[rlregions$source == "dRNH S96"],
  ! is.na(cts$prediction) &
    cts$prediction == "POS" & 
    cts$experiment %in% ip_correct &
    cts$label == "POS" &
    cts$numPeaks > 5000]
ctsmat <- ctsPos@assays@data$cts
vstmat <- DESeq2::vst(ctsmat)
pcdata <- prcomp(t(vstmat))
pcsmall <- pcdata$x[,1:7]
pctpc <- round(100 * pcdata$sdev / sum(pcdata$sdev), digits = 2)
pltdat <- pcsmall %>% 
  as.data.frame() %>%
  rownames_to_column(var = "name") %>%
  inner_join(rlsamples, by = c("name" = "rlsample")) %>%
  mutate(mode = case_when(
    mode %in% RLSeq:::auxdata$mode_cols$mode ~ mode,
    TRUE ~ "misc"
  ))
ip_cols <- c("dRNH" = "#c16fad", "S9.6"="#05bae5")
pltpca <- pltdat %>%
  ggplot(aes(x = PC1, y = PC2, color = ip_type, group = tissue, label = name)) +
  geom_hline(yintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
  geom_point(size=2.5) +
  scale_color_manual(values = ip_cols) +
  theme_bw(base_size = 17)  +
  guides(color=guide_legend(title = NULL))  +
  xlab(paste0("PC1 (", pctpc[1], "%)")) +
  ylab(paste0("PC2 (", pctpc[2], "%)")) 
pltpca
ggsave(plot = pltpca, filename = "results/Figure_5/pca_plot.svg", height = 5, width = 7)

## Get differential RL Regions
cd <- ctsPos@colData %>% 
  as_tibble() %>%
  inner_join(select(rlsamples, rlsample, ip_type), by = c("experiment" = "rlsample")) 
mat <- ctsPos@assays@data$cts
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = mat[,cd$experiment], colData = cd, design = ~ip_type
)
dds <- DESeq2::DESeq(dds)
restbl <- DESeq2::results(dds, contrast = c("ip_type", "S9.6", "dRNH")) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "rlregion") %>%
  as_tibble() %>%
  mutate(rlregion = gsub(rlregion, pattern = "All_", replacement = "")) %>% 
  inner_join(rlregions, by = "rlregion") %>%
  arrange(padj)

## Volcano plot
volcplt <- restbl %>%
  EnhancedVolcano::EnhancedVolcano(
    lab = gsub(.$rlregion, pattern = "All_", replacement = ""), x = "log2FoldChange", y = "padj", 
    pCutoff = 1E-30, FCcutoff = 1,labSize = 4, pointSize = 2,
    title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none"
  ) +
  theme_bw(base_size = 17, base_line_size = 0) +
  ggpubr::rremove("legend")
volcplt
ggsave(plot = volcplt, filename = "results/Figure_5/volcano.png", height = 5.5, width = 7.5)

# Gene set enrichment
ol <- restbl %>%
  mutate(group = ifelse(log2FoldChange > 0, "Over-abundant (S9.6 vs dRNH)", 
                        ifelse(log2FoldChange < 0, "Under-abundant (S9.6 vs dRNH)", "Nothing"))) %>%
  filter(padj < .05) %>%
  group_by(group) %>%
  slice_min(order_by = padj, n = 500) %>%
  mutate(genes = strsplit(allGenes, ",")) %>%
  unnest(genes) %>%
  filter(! is.na(genes)) %>%
  {setNames(group_split(.), nm = group_keys(.)[[1]])} %>%
  lapply(function(x) unique(pull(x, genes)))

resS96 <- enrichR::enrichr(ol$`Over-abundant (S9.6 vs dRNH)`, databases = "GO_Biological_Process_2021")
resdRNH <- enrichR::enrichr(ol$`Under-abundant (S9.6 vs dRNH)`, databases = "GO_Biological_Process_2021")
rescomb <- resS96[[1]] %>%
  as_tibble() %>%
  select(Term, S96__combined_score = Combined.Score, S96__padj = Adjusted.P.value) %>%
  full_join(
    (
      resdRNH[[1]] %>%
        as_tibble() %>% 
        select(Term, dRNH__combined_score = Combined.Score, dRNH__padj = Adjusted.P.value)
    )
  )
pltdata <- rescomb %>% 
  pivot_longer(
    cols = contains("__")
  ) %>%
  mutate(ip_type = gsub(name, pattern = "(.+)__(.+)", replacement = "\\1"),
         metric = gsub(name, pattern = "(.+)__(.+)", replacement = "\\2")) %>%
  select(-name) %>%
  pivot_wider(id_cols = c(Term, ip_type), 
              names_from = metric,
              values_from = value) %>%
  group_by(ip_type) 
terms <- pltdata %>%
  slice_max(order_by = combined_score, n = 6) %>%
  pull(Term)

dd <- pltdata %>%
  filter(Term %in% terms) %>%
  mutate(combined_score = ifelse(! is.na(combined_score), ifelse(combined_score < 1, 0, combined_score), 0)) %>%
  arrange(desc(ip_type), (combined_score)) %>% 
  mutate(
    Term = gsub(Term, pattern = " \\(GO.+", replacement = ""),
    Term = stringr::str_wrap(Term, width = 45)
  )
barplt <- dd %>%
  mutate(
    Term = factor(Term, levels = dd$Term[dd$combined_score > 35]),
    ip_type = gsub(ip_type, pattern = "S96", replacement = "S9.6")
  ) %>%
  ggplot(aes(x = Term, y = combined_score, fill = ip_type)) +
  geom_col(position = position_dodge(.9)) +
  coord_flip() +
  theme_bw(base_size = 15) +
  ylab("Combined Score") +
  xlab(NULL) +
  scale_fill_manual(values = ip_cols) +
  ggpubr::rremove("legend")
barplt
ggsave(barplt, filename = "results/Figure_5/barplt.png", height = 6, width = 10)

# Feature plots
to_plt <- c(
  restbl$rlregion[restbl$log2FoldChange > 0][1:2],
  restbl$rlregion[restbl$log2FoldChange < 0][1:2]
)

plts <- lapply(to_plt, function(pt) {
  pt <- paste0("All_", pt)
  plt <- tibble(
    "Abundance" = vstmat[pt,],
    "name" = names(vstmat[pt,] )
  ) %>% inner_join(pltdat) %>%
    ggplot(aes(x = PC1, y = PC2, color = Abundance, group = tissue)) +
    geom_hline(yintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
    geom_vline(xintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
    geom_point(size=2.5) +
    scale_color_viridis_c(option = "A", direction = -1) +
    theme_bw(base_size = 17) +
    xlab(paste0("PC1 (", pctpc[1], "%)")) +
    ylab(paste0("PC2 (", pctpc[2], "%)")) 
  ggsave(
    plot = plt,
    filename = paste0("results/Figure_5/pca_feature_plot_", gsub(pt, pattern = "All_", replacement = ""),".png"),
    height = 5.5, width = 8
  )
  plt
})
names(plts) <- to_plt

## Explore ribosome results (S5)
gobp <- msigdbr::msigdbr(category = "C2", subcategory = "CP:KEGG")
ribogenes <- gobp$gene_symbol[which(gobp$gs_name == "KEGG_RIBOSOME")]
pltdats <- consInt %>%
  left_join(unique(oltxres), by = c("name", "source")) %>%
  group_by(name) %>%
  mutate(is_ribo = any( SYMBOL %in% ribogenes),
         SYMBOL = paste0(SYMBOL, collapse = "\n")) %>%
  distinct(name, .keep_all = TRUE) %>%
  group_by(source) %>%
  arrange(pct_cons) %>%
  mutate(rank = seq(pct_cons)) %>%
  arrange(desc(pct_cons)) %>%
  mutate(
    label = ifelse(row_number() <= 50 & ! is.na(SYMBOL), SYMBOL, ""),
    label2 = ifelse(row_number() <= 2000 & ! is.na(SYMBOL) & is_ribo, SYMBOL, "")
  )
# riboplt
pltdats %>% 
  arrange(is_ribo) %>%
  ggplot(aes(x = rank, y = pct_cons, color = is_ribo, label=label2)) +
  geom_hline(yintercept = 50, color = "grey", linetype = "dashed", alpha = .85) +
  geom_point() +
  geom_rug(sides = "b") +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "#e6e6e6")) +
  # ggrepel::geom_label_repel(max.overlaps = Inf, color = "black", 
  # force = 50, nudge_x = -500, nudge_y = -.05) +
  facet_wrap(~ source, scales = "free_x") +
  theme_bw(base_size = 14) +
  ylab("RL region conservation (% samples)") +
  xlab("RL region rank") +
  guides(color = guide_legend(title = "Ribosomal Protein"))

# Get venn
toplt <- pltdats %>%
  slice_max(n = 2000, order_by = pct_cons) %>%
  mutate(SYMBOL = strsplit(SYMBOL, split = "\n")) %>%
  unnest(SYMBOL) %>%
  filter(SYMBOL %in% ribogenes)
complst <- data.frame(
  "S9.6" = ribogenes %in% unique(toplt$SYMBOL[toplt$source == "S9.6"]),
  "dRNH" = ribogenes %in% unique(toplt$SYMBOL[toplt$source == "dRNH"]),
  "Ribo genes" = rep(TRUE, length(ribogenes))
) %>% as.matrix()
fit <- eulerr::euler(complst)
pdf("results/Figure_5/euler_ribo.pdf")
plot(fit,
     quantities = TRUE,
     lty = 1, fills = c("#0dbae5", "#c170ad", "white"),
     labels = list(font = 16))
dev.off()

# Pick which gene to show on genome browser (top diff where dRNH does not equal 0)
pltdats %>% 
  mutate(SYMBOL = strsplit(SYMBOL, split = "\n")) %>%
  unnest(SYMBOL) %>%
  filter(SYMBOL %in% ribogenes) %>%
  select(SYMBOL, pct_cons, source) %>%
  group_by(SYMBOL, source) %>%
  summarize(pct_cons = max(pct_cons)) %>%
  pivot_wider(id_cols = SYMBOL, values_from = pct_cons, names_from = source, values_fill = 0) %>%
  mutate(diff = S9.6 - dRNH) %>%
  arrange(desc(diff))


##### Figure 6: dRNH-only mapping of enhancers #####

resloc <- "results/Figure_6/"
dir.create("results/Figure_9", showWarnings = FALSE)

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
      "dRNH-shared" = "#ccb6cc", 
      "S9.6" = "#82d0e8"
    )
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("dRNH-only", "dRNH-shared"))
  ) +
  ggpubr::rremove("legend")
plt
ggsave(plt, filename = file.path(resloc, "consensus_peak_width.png"))


### Proportion test
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
    group3 = factor(group3, levels = (unique(group3))),
    group2 = factor(group2, levels = rev(unique(group2)))
  ) %>%
  group_by(group2) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = group2, y = prop, fill = group3)) +
  geom_col(color="black") +
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  scale_fill_manual(
    values = c(
      "Unshared" = "#c4c4c4",
      "Shared" = "#d88a19"
    )
  ) +
  ggtitle("Consensus sites shared vs unshared")
plt
ggsave(plt, filename = file.path("results/Figure_4/", "prop_consensus_shared_barchart.png"), height = 6, width = 8)

### Analysis of Tx Features
## What are the features overlapping consensus peaks in each group?
## How does this relate to peak size/strength?
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
num_sel <- 6
db <- "KEGG_2021_Human"
db <- "ChEA_2016"
db <- "MSigDB_Hallmark_2020"
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

### cCREs
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
    angle_col = 45, fontsize = 17,
    main = "Enrichment of peaks within CREs",
    cluster_cols = F
  )
png(file.path(resloc, "CCRE_heatmap.png"))
plt
dev.off()

#### Analysis of dRNH-only peaks and enhancers 

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

## Step 0: Get all the enhancers
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
  ChIPseeker::annotatePeak(TxDb = txdb, tssRegion = c(-3000, 3000), verbose = T)
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
  geom_col(color="black") +
  coord_flip() +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  scale_fill_manual(
    values = c(
      "non-Enhancer" = "#d1d1d1",
      "Enhancer" = "#2b41af"
    ), labels=c("Non-enhancer", "Enhancer")
  ) +
  guides(fill=guide_legend(title = NULL)) +
  ggtitle("Consensus peak proportion at enhancers") +
  theme_bw(base_size = 16)
plt
ggsave(plt, filename = file.path(resloc, "consensus_pks_prop_all_enhancers.png"))

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
  geom_col(color="black") +
  coord_flip() +
  theme_gray(base_size = 14) +
  theme(legend.title = element_blank()) +
  ylab("Proportion of consensus peaks") +
  xlab(NULL) +
  scale_fill_manual(
    values = c(
      "non-distal-enhancer" = "#d1d1d1",
      "distal-enhancer" = "#2b41af"
    ), labels=c("Other", "Distal Enhancer")
  ) +
  guides(fill=guide_legend(title = NULL)) +
  ggtitle("Consensus peak proportion at distal enhancers") +
  theme_bw(base_size = 16)

plt
ggsave(plt, filename = file.path(resloc, "consensus_pks_denh.png"))

##### R2R: dEnh R-loop size #####
pkenhsize <- lapply(pks, function(x) {
  y <- gr2bed2(x)
  ints <- valr::bed_intersect(y, gr2bed2(dEnh))
  nms <- unique(y$name)
  tibble(
    "names" = nms
  ) %>%
    mutate(inenh = ifelse(names %in% ints$name.x, "Distal Enhancer", "Other")) %>%
    group_by(inenh)
}) %>% bind_rows(.id = "group")

data4plot <- cons2 %>%
  bind_rows(.id = "group") %>% 
  inner_join(pkenhsize, by = c("group", "name" = "names")) %>% 
  mutate(width = end - start) %>% 
  select(group, name, width, inenh)


plt <- data4plot %>% 
  ggplot(aes(x = inenh, y = width, fill = inenh)) +
  geom_violin(width = .6, position = position_dodge(.5), trim = F, alpha = .3) +
  geom_boxplot(width = .2, position = position_dodge(.5), alpha = .3) +
  geom_jitter(position = position_jitterdodge(.1), size = .25, alpha = .2) +
  facet_wrap(~group) + 
  scale_y_log10(limits = c(1E2, 5E5)) +
  theme_bw(base_size = 20) +
  ylab("Peak width (log bp)") +
  xlab(NULL) +
  ggtitle("Consensus peak width distribution") +
  guides(
    fill = guide_legend(title = NULL), 
    color = guide_legend(title = NULL)
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Distal Enhancer", "Other")),
    label = "p.format",
    size = 5,
  ) +
  scale_fill_manual(
    values = c(
      "Other" = "#d2d2d2",
      "Distal Enhancer" = "#38499f"
    ), labels=c("Other", "Distal Enhancer")
  ) +
  ggpubr::rremove("legend") 
plt

#################################


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
peakAnnoList <- lapply(
  pks,
  ChIPseeker::annotatePeak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000), 
  verbose = T
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

##### Figure 7: iPSC/CUTLL1 distal enhancers #####

dat <- RLHub::rlbase_samples()
resloc <- "results/Figure_7"
dir.create(resloc, showWarnings = FALSE)

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
  grs, ChIPseeker::annotatePeak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000), verbose = T
)
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)
gri <- reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]

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
num_sel <- 6
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
  theme_bw(base_size = 20) +
  ylab("Combined Score") +
  xlab(NULL) +
  ggtitle("CUTLL1 dRNH-enhancer gene targets", subtitle = "Enrichment in CellMarker Database")
plt
ggsave(plt, filename = file.path(resloc, "Cellmarker_CUTLL1_dENH.png"))


### CTCF
CTCF_CUTLL1_1 <- valr::read_narrowpeak(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3732nnn/GSM3732747/suppl/GSM3732747_CUTLL1_DMSO_sort_peaks.narrowPeak.bed.gz"
) 
CTCF_CUTLL1_2 <- valr::read_narrowpeak(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3732nnn/GSM3732744/suppl/GSM3732744_CUTLL1_CTCF2_sort_peaks.narrowPeak.bed.gz"
)
CTCF_CUTLL1_3 <- valr::read_narrowpeak(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3732nnn/GSM3732743/suppl/GSM3732743_CUTLL1_CTCF1_sort_peaks.narrowPeak.bed.gz"
)
CTCF_CUTLL1 <- reduce(
  c(
    GenomicRanges::makeGRangesFromDataFrame(CTCF_CUTLL1_1),
    GenomicRanges::makeGRangesFromDataFrame(CTCF_CUTLL1_2),
    GenomicRanges::makeGRangesFromDataFrame(CTCF_CUTLL1_3)
  )
)

pltd <- ChIPpeakAnno::binOverFeature(
  gris[[1]], gris[[2]], gris[[3]], gris[[4]], gris[[5]], gris[[6]],
  annotationData = CTCF_CUTLL1,
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
  ggtitle("Peak pileup around CTCF peaks", subtitle = "CUTLL1 MapR (dRNH) peaks") +
  ylab("Peak density") +
  xlab("Distance to CTCF (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2")
plt

# Get the intergenic peaks from both CUTTL1 samples & find overlap
gri <- reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]
olde32 <- ChIPpeakAnno::findOverlapsOfPeaks(
  GenomicRanges::makeGRangesFromDataFrame(CTCF_CUTLL1) %>% reduce(),
  gri
)
# PCT olap
length(unique(olde32$overlappingPeaks$`GenomicRanges..makeGRangesFromDataFrame.CTCF_CUTLL1......reduce..///gri`$peaks2)) / length(unique(olde31$all.peaks$gri))
ChIPpeakAnno::makeVennDiagram(
  olde32,
  NameOfPeaks = c("CTCF peaks", "In"),
  fill=c("gray", "#0BB490"), margin=.1
)


### SMC3
SMC3_CUTLL1_1 <- valr::read_narrowpeak(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4230nnn/GSM4230108/suppl/GSM4230108_SMC3_rep1_sort_peaks.narrowPeak.bed.gz"
) 
SMC3_CUTLL1_2 <- valr::read_narrowpeak(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4230nnn/GSM4230109/suppl/GSM4230109_SMC3_rep2_sort_peaks.narrowPeak.bed.gz"
)
SMC3_CUTLL1 <- reduce(
  c(
    GenomicRanges::makeGRangesFromDataFrame(SMC3_CUTLL1_1),
    GenomicRanges::makeGRangesFromDataFrame(SMC3_CUTLL1_2)
  )
)

pltd <- ChIPpeakAnno::binOverFeature(
  gris[[1]], gris[[2]], gris[[3]], gris[[4]], gris[[5]], gris[[6]],
  annotationData = SMC3_CUTLL1,
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
  ggtitle("Peak pileup around SMC3 peaks", subtitle = "CUTLL1 MapR (dRNH) peaks") +
  ylab("Peak density") +
  xlab("Distance to SMC3 (bp) (5->3)") +
  theme_bw(14) +
  theme(legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2")
plt

# Get the intergenic peaks from both iPSC samples & find overlap
gri <- reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]
olde33 <- ChIPpeakAnno::findOverlapsOfPeaks(
  GenomicRanges::makeGRangesFromDataFrame(CTCF_CUTLL1) %>% reduce(),
  GenomicRanges::makeGRangesFromDataFrame(SMC3_CUTLL1) %>% reduce(),
  gri
)
# PCT olap
length(unique(olde31$overlappingPeaks$`GenomicRanges..makeGRangesFromDataFrame.CTCF_iPS......reduce..///gri`$peaks2)) / length(unique(olde31$all.peaks$gri))
ChIPpeakAnno::makeVennDiagram(
  olde33,
  NameOfPeaks = c("CTCF peaks", "SMC3 peaks", "In"),
  fill=c("gray", "blue", "#0BB490"), margin=.1
)
combo <- c(A = 23189, B = 110, C = 3320, "A&B" = 4425, "A&C" = 963, "B&C" = 7, "A&B&C" = 289)
fit <- eulerr::euler(
  combo  
)
plot(fit)

# Proportion plot
tibble(
  group = c(
    "SMC3 + CTCF", "CTCF-only", "SMC3 + CTCF", "CTCF-only"
  ),
  group2 = c("No MapR", "No MapR", "MapR", "MapR"),
  size = c(
    4425, 23189, 289, 963
  )
) %>% 
  group_by(group2) %>% 
  mutate(pct = size / sum(size)) %>% 
  ggplot(aes(x = group2, y = pct, fill = group)) +
  geom_bar(stat = "identity", position = "stack", color="black") +
  theme_bw(base_size = 14) +
  coord_flip() +
  ylab("Peak colocalization (%)") +
  xlab(NULL) +
  scale_fill_manual(
    values = c(
      "CTCF-only" = "#D1D1D2",
      "SMC3 + CTCF" = "#9F87B5"
    ),
    guide = guide_legend(title = NULL, reverse = TRUE)
  ) 





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
  grs, ChIPseeker::annotatePeak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000), verbose = T
)

## Find overlaps with enhancers

# Get intergenic ranges
gris <- lapply(seq(pal), function(i) {
  grs[[i]][pal[[i]]@detailGenomicAnnotation$Intergenic, ]
})
names(gris) <- names(pal)
gri <- reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]

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
  theme_bw(20) +
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
  theme_bw(base_size = 20) +
  ylab("Combined Score") +
  xlab(NULL) +
  ggtitle("iPSCs dRNH-enhancer gene targets", subtitle = "Enrichment in CellMarker Database")
plt


#### Examine SA1/2/CTCF and TAD boundaries

### CTCF
CTCF_iPS <- valr::read_narrowpeak(
  "https://www.encodeproject.org/files/ENCFF322WKG/@@download/ENCFF322WKG.bed.gz"
) %>% 
  filter(qvalue > 4)
pltd <- ChIPpeakAnno::binOverFeature(
  gris[[1]], gris[[2]], gris[[3]], gris[[4]],
  annotationData = GenomicRanges::makeGRangesFromDataFrame(CTCF_iPS),
  radius = 5000, nbins = 100, FUN = length, errFun = 0,
  # featureSite = "bothEnd",
  # PeakLocForDistance = "middle",
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
  ggtitle("Peak pileup around CTCF peaks", subtitle = "iPSCs MapR (dRNH) peaks") +
  ylab("Peak density") +
  xlab("Distance to CTCF (bp) (5->3)") +
  theme_bw(20) +
  theme(legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2")
plt

# Get the intergenic peaks from both iPSC samples & find overlap
gri <- reduce(do.call("c", unlist(gris, use.names = F)), with.revmap = T)
gri <- gri[which(sapply(gri$revmap, length) > 1), ]
olde31 <- ChIPpeakAnno::findOverlapsOfPeaks(
  GenomicRanges::makeGRangesFromDataFrame(CTCF_iPS) %>% reduce(),
  gri
)
# PCT olap
length(unique(olde31$overlappingPeaks$`GenomicRanges..makeGRangesFromDataFrame.CTCF_iPS......reduce..///gri`$peaks2)) / length(unique(olde31$all.peaks$gri))
ChIPpeakAnno::makeVennDiagram(
  olde31,
  NameOfPeaks = c("CTCF peaks", "Intergenic dRNH sites"),
  fill=c("gray", "#F58776"), margin=.1
)




### Loop anchors
a2a <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3177nnn/GSM3177704/suppl/GSM3177704_iPS.loops.txt.gz",
  col_names = c("a1", "a2", "score", "signal", "pval")
)
locpat <- "(.+):(.+)\\-(.+)"
a2a <- a2a %>% 
  mutate(lid = row_number()) %>% 
  pivot_longer(cols = c(a1, a2)) %>% 
  mutate(aid = paste0(name, "_", lid)) %>% 
  mutate(
    chrom = gsub(value, pattern = locpat, replacement = "\\1"),
    start = gsub(value, pattern = locpat, replacement = "\\2") %>% as.numeric(),
    end = gsub(value, pattern = locpat, replacement = "\\3") %>% as.numeric(),
    pval = ifelse(pval == 0, .Machine$double.xmin, pval),
    pval = -log10(pval)
  ) %>% 
  select(chrom, start, end, score=pval, name=aid)
topanchor <- a2a %>% 
  filter(score > 10)

# liftover
chain <- "data/hg19_to_hg38.chain.gz"
if (!file.exists(gsub(chain, pattern = "\\.gz", replacement = ""))) {
  download.file(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz", 
    destfile = chain
  )
  R.utils::gunzip(chain)
}
chn2 <- rtracklayer::import.chain(gsub(chain, pattern = "\\.gz", replacement = ""))
lo <- GenomicRanges::makeGRangesFromDataFrame(topanchor, keep.extra.columns = T) %>% 
  rtracklayer::liftOver(chain = chn2) %>% 
  unlist() %>% 
  unique()

download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5032nnn/GSM5032462/suppl/GSM5032462_KO_D0_iPS_Rep1_RawMatrix_1kb.h5",
  destfile = "tmp/iPS_HiC.h5"
)
iphic <- rhdf5::h5read("tmp/iPS_HiC.h5", name = "intervals") %>%
  bind_cols() %>% 
  select(chrom = chr_list, start = start_list, end = end_list, score = extra_list) %>% 
  filter(score > 5)
iphichg38 <- iphic %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>% 
  rtracklayer::liftOver(chain = chn2) %>% 
  unlist() %>% 
  unique()


### SA1/2 & consensus

## Get CTCF consensus 

annCTCF <- annFull$encodeTFBS__CTCF
pks4gr <- lapply(pks4, GenomicRanges::makeGRangesFromDataFrame)
pltd <- ChIPpeakAnno::binOverFeature(
  pks4gr[[1]], pks4gr[[2]], pks4gr[[3]],
  annotationData = GenomicRanges::makeGRangesFromDataFrame(annCTCF),
  radius = 5000, nbins = 100, FUN = length, errFun = 0,
  # featureSite = "bothEnd",
  # PeakLocForDistance = "middle",
  xlab = "Distance from Enh (bp)", ylab = "count"
)
plt <- pltd %>%
  as.data.frame() %>%
  rownames_to_column("pos") %>%
  as_tibble() %>%
  mutate(pos = as.numeric(pos)) %>%
  dplyr::rename("dRNH-shared" = 3, "dRNH-only" = 2, "S9.6" = 4) %>%
  pivot_longer(cols = !contains("pos")) %>%
  ggplot(aes(x = pos, y = value, color = name)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .25) +
  geom_point(size=1.6) +
  geom_line(size=1.1) +
  ggtitle("Peak pileup around CTCF consensus", subtitle = "RL consensus peaks") +
  ylab("Peak density") +
  xlab("Distance to CTCF consensus peak (bp) (5->3)") +
  theme_bw(20) +
  theme(legend.title = element_blank()) +
  scale_color_manual(
    values = c(
      "dRNH-only" = "#e2a3e3", 
      "dRNH-shared" = "#ccb6cc", 
      "S9.6" = "#82d0e8"
    )
  )
plt

## Full Cohesin

annSA1 <- annFull$Cohesin__STAG1
annSA2 <- annFull$Cohesin__STAG2
annSMC3 <- annFull$encodeTFBS__SMC3
annRAD21 <- annFull$encodeTFBS__RAD21
list(
  "SA1" = annSA1,
  "SA2" = annSA2,
  "RAD21" = annRAD21,
  "SMC3" = annSMC3,
  "CTCF" = annCTCF
) -> ylst
gen <- valr::read_genome("tmp/hg38.chrom.sizes")

plt <- pks4 %>% lapply(
  function(x) {
    x <- valr::gr_to_bed(x)
    lapply(
      ylst, function(y) {
        set.seed(42); valr::bed_fisher(
          x,
          # Down sample to match size range of other annotations
          y = slice_sample(y, n = 2e4), 
          genome = gen
        )
      }
    ) %>% bind_rows(.id = "group")
  }
) %>% 
  bind_rows(.id = "group2") %>% 
  mutate(
    p.value = ifelse(p.value == 0, .Machine$double.xmin, p.value),
    `pval (-log10)` = -log10(p.adjust(p.value))
  ) %>% 
  dplyr::rename(`Odds Ratio`=estimate) %>% 
  mutate(
    group2 = factor(group2, levels = rev(unique(.$group2)))
  ) %>% 
  ggplot(aes(x = group2, y = `pval (-log10)`, fill = `Odds Ratio`)) +
  geom_col(color = "black") +
  facet_wrap(~group, ncol = 5) +
  theme_bw(16) +
  ggtitle("Peak enrichment at CTCF/Cohesin", subtitle = "RL consenus peaks") +
  coord_flip() +
  xlab(NULL) +
  ylab("Fisher test padj (-log10)") 
plt

########################



#### eRNA analysis ####

fls <- list(
  "HEK293" = "https://www.encodeproject.org/files/ENCFF476TTU/@@download/ENCFF476TTU.bed.gz",
  "HCT116" = "https://www.encodeproject.org/files/ENCFF513PJK/@@download/ENCFF513PJK.bed.gz",
  "iPSCs" = "https://www.encodeproject.org/files/ENCFF115RIR/@@download/ENCFF115RIR.bed.gz"
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

chain <- "data/hg38_to_hg19.chain.gz"
if (!file.exists(gsub(chain, pattern = "\\.gz", replacement = ""))) {
  download.file(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz", 
    destfile = chain
  )
  R.utils::gunzip(chain)
}
chn <- rtracklayer::import.chain(gsub(chain, pattern = "\\.gz", replacement = ""))

### iPSCs

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

## Get dEnh2
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

# Get intersect of GRO
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
    label = "p.format", size = 5
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_classic(base_size = 20) +
  scale_fill_manual(
    values = c(
      "dRNH-accessible dEnh" = "#F58776",
      "Total dEnh" = "grey"
    )
  ) +
  ggpubr::rremove("legend")

plt
ggsave(plt, filename = file.path(resloc, "iPSC_eRNA_expression_dENH.png"))



### CUTLL1
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

gro <- valr::bed_merge(
  bind_rows(list(
    gro1,
    gro2,
    gro3,
    gro4
  )
  ), score = sum(score))

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
    label = "p.format", size = 5
  ) +
  xlab(NULL) +
  ylab("eRNA RPKM (log2 + 1)") +
  ggtitle(
    paste0(ct, " distal enhancer RNA profile"),
    subtitle = "dRNH-accessible vs total enhancer population"
  ) +
  theme_classic(base_size = 20) +
  scale_fill_manual(
    values = c(
      "dRNH-bound dEnh" = "#0BB490",
      "Total dEnh" = "gray"
    )
  ) +
  ggpubr::rremove("legend")
plt
ggsave(plt, filename = file.path(resloc, "CUTLL1_eRNA_expression_dENH.png"))


##### Figure 8: Bivalent enhancers iPSCs #####

resloc <- "results/Figure_8"
dir.create(resloc, showWarnings = FALSE)
num_sel <- 6
chain <- "data/hg38_to_hg19.chain.gz"
chn <- rtracklayer::import.chain(gsub(chain, pattern = "\\.gz", replacement = ""))

## hg19 windows
download.file("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes", destfile = "tmp/hg19.chrom.sizes")
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
  theme_bw(base_size = 16) +
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

lenmin <- 1  # Clusters must contain more than lenmin enhancers

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


#### Analysis of eCLIP/ChIP binding

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

### Comparison plot
rbind(mutate(dd3, group = "Genic"), mutate(dd2, group = "Intergenic")) %>% 
  pivot_longer(cols = contains("-")) %>% 
  mutate(name2 = str_c(name, group, sep = " - ")) %>% 
  group_by(group, name, ip) %>% 
  filter(! is.na(value)) %>% 
  summarise(mn=mean(value)) %>% 
  filter(ip == "eCLiP") %>% 
  pivot_wider(id_cols = name, names_from = group, values_from = mn) %>% 
  mutate(dif = Genic / Intergenic)
plt <- rbind(mutate(dd3, group = "Genic"), mutate(dd2, group = "Intergenic")) %>% 
  pivot_longer(cols = contains("-")) %>% 
  mutate(name2 = str_c(name, group, sep = " - ")) %>% 
  # filter(ip == "ChIP") %>% 
  ggplot(aes(x = name2, y = value, color = group)) +
  geom_violin(trim = F) +
  geom_jitter(width = .15) +
  facet_wrap(~ip) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("S9.6-only - Genic", "S9.6-only - Intergenic"),
      c("dRNH-only - Genic", "dRNH-only - Intergenic")
    )
  ) +
  theme_bw(14) +
  scale_y_continuous(limits=c(-.5, 8), expand = c(0,0)) +
  xlab(NULL) +
  ggpubr::rotate_x_text(30) +
  ylab("Enrichment (log2 odds ratio)") +
  ggtitle("Enrichment of consensus peaks within ChIP/eCLiP sites")
plt
ggsave(plt, filename = file.path(resloc, "eCLiP_vs_ChIP_binding.png"))


##### Figure 9: Conserved vs Variable #####
dir.create("results/Figure_9", showWarnings = FALSE)

## Define conserved vs variable
totalrnh <- length(unique(rlsamples$rlsample[rlsamples$ip_type == "dRNH" & rlsamples$label == "POS" & rlsamples$genome == "hg38" & rlsamples$numPeaks != 0]))
totalS96 <- length(unique(rlsamples$rlsample[rlsamples$ip_type == "S9.6" & rlsamples$label == "POS" & rlsamples$genome == "hg38" & rlsamples$numPeaks != 0]))
# Apply scaling factor
rlregions_cons <- rlregions %>%
  dplyr::rename(
    pct_cons = conservation_pct
  ) %>%
  arrange(desc(pct_cons)) %>%
  mutate(rank = row_number()) %>%
  mutate(quant = case_when(
    pct_cons > 60 ~ "60-100%",
    pct_cons > 40 ~ "40-60%",
    pct_cons > 26 ~ "26-40%",
    pct_cons > 20 ~ "20-26%",
    TRUE ~ "15-20%"
  ))

## Summit the peaks (we need to transfer summits from the dRNH and S9.6 consensus sites)
consSum <- consInt %>%
  select(chrom, so, eo, source, start, end, peak, name) %>%
  distinct(name, .keep_all = T)
locpat <- "(.+):(.+)-(.+):(.+)"
rlregions_cons <- rlregions_cons %>%
  mutate(
    chrom = gsub(location, pattern = locpat, replacement = "\\1"),
    so = as.numeric(gsub(location, pattern = locpat, replacement = "\\2")),
    eo = as.numeric(gsub(location, pattern = locpat, replacement = "\\3")),
    strand = gsub(location, pattern = locpat, replacement = "\\4")
  )
rljoin <- inner_join(rlregions_cons, consSum, by = c("chrom", "so")) %>%
  mutate(source = source.x)
rljoin <- rljoin %>% distinct(rlregion, .keep_all = TRUE)

## Percent conservation rank plot
quant_cols <- c(
  "60-100%" = "#5e03fc",
  "40-60%" = "#a678f5",
  "26-40%" = "#a89ff5",
  "20-26%" = "#bdcdff",
  "15-20%" = "#c3e4f7"
)
plt <- rljoin %>%
  ggplot(aes(x = rank, y = pct_cons, color = quant)) +
  geom_hline(yintercept = 50, color = "grey", linetype = "dashed") +
  geom_point() +
  scale_x_reverse() +
  geom_rug(
    # sides = "b"
  ) +
  ylab("Conservation across samples (%)") +
  xlab("RL region rank") +
  scale_color_manual(values = quant_cols) +
  guides(color = guide_legend(title = "Percentile range")) +
  theme_bw(base_size = 16)
ggsave(plt, filename = "results/Figure_9/pct_cons_ranks.svg", height = 6, width = 9)

## Ridge plot of conservation
plt <- rlregions_cons %>%
  arrange(desc(pct_cons)) %>%
  mutate(source = ifelse(source == "dRNH S96", "Both", ifelse(source == "S96", "S9.6", "dRNH")),
         source = factor(source, levels = c("S9.6", "dRNH", "Both"))) %>%
  ggplot(aes(x = pct_cons, y = source, color = source, fill = source)) +
  ggridges::geom_density_ridges(alpha = .6) +
  ggridges::theme_ridges(font_size = 16) +
  scale_x_log10(breaks = c(15, 20, 30, 50, 100)) +
  scale_fill_manual(
    values = c("Both" = "#6795C9", "dRNH" = "#C170AD", "S9.6" = "#0DBAE5")
  ) +
  scale_color_manual(
    values = c("Both" = "#6795C9", "dRNH" = "#C170AD", "S9.6" = "#0DBAE5")
  ) +
  xlab("Conservation (%, log10)") +
  ylab(NULL)
ggsave(plt, filename = "results/Figure_9/pct_cons_ridges.svg")

## Get the overlap stats with genomic
oltx <- rljoin %>%
  select(chrom, start, end, rlregion, source) %>%
  valr::bed_intersect(txfeats)
oltxsum <- oltx %>%
  group_by(rlregion.x) %>%
  summarise(
    type = ifelse("TSS" %in% type.y, "TSS", 
                  ifelse("TTS" %in% type.y, "TTS",
                         ifelse("fiveUTR" %in% type.y, "fiveUTR", 
                                ifelse("threeUTR" %in% type.y, "threeUTR", 
                                       ifelse("Exon" %in% type.y, "Exon", "Intron")))))
  ) %>%
  dplyr::rename(rlregion = rlregion.x)
oltxres <- oltxsum %>%
  full_join(rljoin, by = "rlregion") %>%
  mutate(type = ifelse(is.na(type), "Intergenic", type)) %>%
  distinct(rlregion, .keep_all = TRUE)
plt <- oltxres %>%
  group_by(quant, type) %>%
  tally() %>%
  mutate(n = 100*n / sum(n),
         type = factor(type, levels = rev(c(
           "TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "Intergenic"
         ))),
         quant = factor(quant, levels = rev(unique(quant))))  %>%
  ggplot(aes(x = quant, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  coord_flip() +
  ylab("Peak location (%)") +
  xlab(NULL) +
  scale_fill_manual(values = txfeat_cols,
                    guide = guide_legend(title = "Feature", reverse = TRUE)) 
ggsave(plt, filename = "results/Figure_9/pct_cons_txfeats.svg")

## Split regions by quantiles of conservation
## Do Enrichment analysis
# tssrls <- oltxres$rlregion[! oltxres$type %in% c("Intron", "Intergenic")]
rllst <- rljoin %>%
  group_by(quant) %>%
  {setNames(group_split(.), nm = group_keys(.)[[1]])}

## Make barplots
rllstenr <- pbapply::pblapply(
  rllst,
  function(x) {
    gen <- x %>%
      # filter(rlregion %in% tssrls) %>%
      pull(allGenes)
    genlst <- unique(unlist(strsplit(gen, ",")))
    enrichR::enrichr(genes = genlst, databases = c("GO_Biological_Process_2021", "ChEA_2016",
                                                   "KEGG_2021_Human", "MSigDB_Hallmark_2020"))
  }
)

# Check XRN2 termination
xrn2genes <- unlist(strsplit(rllstenr$`60-100%`$ChEA_2016$Genes[1], ";"))
rllst$`60-100%` %>%
  mutate(allGenes = strsplit(allGenes, ",")) %>%
  unnest(cols = allGenes) %>%
  filter(allGenes %in% xrn2genes)

num_sel <- 5
db <- "KEGG_2021_Human"
db <- "ChEA_2016"
db <- "MSigDB_Hallmark_2020"
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
      combined_score > 80 ~ 80,
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
  ylab("RL Region conservation percentile") +
  xlab(NULL) +
  scale_color_viridis_c(direction = -1, option = "A", end = .9,
                        guide = guide_colorbar(title = "P adj (-log10)")) 
ggsave(plt, filename = "results/Figure_9/path_enrich_tf.svg", height = 8, width = 6)

## Housekeeping genes
# msigdbr_df <- msigdbr::msigdbr(category = "C2", subcategory = "CGP")
# msigdbr_t2g <- msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
hskg <- msigdbr::msigdbr(category = "C2", subcategory = "CGP") %>%
  filter(gs_name == "HSIAO_HOUSEKEEPING_GENES") %>% pull(gene_symbol)
allgen <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
  AnnotationDbi::select(
    ., keys = AnnotationDbi::keys(.), columns = "SYMBOL"
  ) %>% pull(SYMBOL)
rllstgenes <-sapply(
  rllst,
  function(x) {
    gen <- x %>%
      # filter(rlregion %in% tssrls) %>%
      pull(allGenes)
    genlst <- unique(unlist(strsplit(gen, ",")))
    print(calculate.overlap.and.pvalue(genlst, hskg, lower.tail = FALSE, total.size = length(allgen)))[3]
  }
)
plt <- tibble(
  "cons_pct" = names(rllstgenes),
  "enrichment" = rllstgenes
) %>%
  mutate(enrichment = -log10(p.adjust(enrichment))) %>%
  ggplot(aes(x = cons_pct, y  = enrichment)) +
  geom_col(fill = "#57b580") +
  coord_flip() + ylab("Enrichment [-log10(padj)]") + 
  xlab("RL region conservation percentile") +
  theme_bw(base_size = 14)
ggsave(
  plt,
  filename = "results/Figure_9/housekeeping_genes.svg",
  height = 5, width = 5
)

## Table S2 -- RL regions
dir.create("results/Table_S2/", showWarnings = FALSE)
locpat <- "(.+):(.+)-(.+):(.+)"
rljoin %>%
  mutate(
    chrom = gsub(location, pattern = locpat, replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = locpat, replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = locpat, replacement = "\\3")),
    rlregion = gsub(rlregion, pattern = "All_", replacement = "")
  ) %>%
  select(
    RLRegion = rlregion,
    chrom, start, end, 
    `Source (ip_type)` = source,
    `Conservation (%)` = pct_cons,
    `# Studies` = nStudies,
    `Genes (overlapping)` = allGenes,
    `In repeat region` = is_repeat,
    `In RLFS` = is_rlfs
  ) %>% 
  write_csv("results/Table_S2/rlregions.csv")

##### Table S3 Differentially-abundant RL Regions
dir.create("results/Table_S3/", showWarnings = FALSE)
restbl %>%
  mutate(rlregion = gsub(rlregion, pattern = "All_", replacement = "")) %>%
  select(rlregion, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  write_csv("results/Table_S3/differentially_abundant_rlregions.csv")




### Pausing index


## HeLa
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5420nnn/GSM5420820/suppl/GSM5420820_Async_repA_PROseq.bed.gz"
df <- gsub(x, pattern = ".+/suppl/(.+)", replacement = "tmp/\\1")
if (! file.exists(df)) download.file(x, destfile = df)
## POS
infile <- gsub(df, pattern = "\\.bed\\.gz", replacement = ".pos.bed") 
infile2 <- gsub(df, pattern = "\\.bed\\.gz", replacement = ".pos.sort.bed") 
outfile <-  gsub(df, pattern = "\\.bed\\.gz", replacement = ".pos.bdg")
outfile2 <-  gsub(df, pattern = "\\.bed\\.gz", replacement = ".pos.bw") 
if (! file.exists(outfile2)) {
  y <- valr::read_bed12(x)
  chroms <- valr::read_genome("tmp/hg38.chrom.sizes")
  y %>% 
    filter(strand == "+", chrom %in% chroms$chrom) %>%
    write_tsv(file = infile, col_names = FALSE)
  system(paste0("tail -n +2 ", infile, " | sort -k1,1 > ", infile2))
  system(paste0("~/miniconda3/condabin/conda run -n rlbaseData genomeCoverageBed -bg -i ", infile2, " -g tmp/hg38.chrom.sizes > ", outfile))
  system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedGraphToBigWig ", outfile, " tmp/hg38.chrom.sizes ", outfile2))
}
aws.s3::put_object(file = outfile2, object = paste0("misc/", outfile2), bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)

## NEG
infile <- gsub(df, pattern = "\\.bed\\.gz", replacement = ".neg.bed") 
infile2 <- gsub(df, pattern = "\\.bed\\.gz", replacement = ".neg.sort.bed") 
outfile <-  gsub(df, pattern = "\\.bed\\.gz", replacement = ".neg.bdg")
outfile2 <-  gsub(df, pattern = "\\.bed\\.gz", replacement = ".neg.bw") 
if (! file.exists(outfile2)) {
  y %>% 
    filter(strand == "-", chrom %in% chroms$chrom) %>%
    write_tsv(file = infile, col_names = FALSE)
  system(paste0("tail -n +2 ", infile, " | sort -k1,1 > ", infile2))
  system(paste0("~/miniconda3/condabin/conda run -n rlbaseData genomeCoverageBed -bg -i ", infile2, " -g tmp/hg38.chrom.sizes > ", outfile))
  system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedGraphToBigWig ", outfile, " tmp/hg38.chrom.sizes ", outfile2))
}
aws.s3::put_object(file = outfile2, object = paste0("misc/", outfile2), bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)



## B2B
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5103nnn/GSM5103029/suppl/GSM5103029_B2BPROseq_V1.bedGraph.gz"
df <- gsub(x, pattern = ".+/suppl/", replacement = "")
if (! file.exists(df)) download.file(x, destfile = df, method = "curl")
y <- read_tsv(df, col_names = c("chrom", "start", "end", "score"))
## POS
infile <- gsub(df, pattern = "\\.bedGraph\\.gz", replacement = ".pos.bdg") %>% file.path("tmp", .)
outfile <- gsub(df, pattern = "\\.bedGraph\\.gz", replacement = ".pos.bw") %>% file.path("tmp", .)
chroms <- valr::read_genome("tmp/hg38.chrom.sizes")
y %>% filter(
  score > 0,  chrom %in% chroms$chrom
) %>% 
  write_tsv(file = infile, col_names = FALSE)
system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedGraphToBigWig ", infile, " tmp/hg38.chrom.sizes ", outfile))
## NEG
infile <- gsub(df, pattern = "\\.bedGraph\\.gz", replacement = ".neg.bdg") %>% file.path("tmp", .)
outfile <- gsub(df, pattern = "\\.bedGraph\\.gz", replacement = ".neg.bw") %>% file.path("tmp", .)
y %>% filter(
  score < 0,   chrom %in% chroms$chrom
) %>% 
  mutate(score = -1 * score) %>% 
  write_tsv(file = infile, col_names = FALSE)
system(paste0("~/miniconda3/condabin/conda run -n rlbaseData bedGraphToBigWig ", infile, " tmp/hg38.chrom.sizes ", outfile))


### Get HCT116
download.file("https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz", destfile = "tmp/hg19_hg38.chain.gz")
x <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156737/suppl/GSE156737_WT_DMSO.ps.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))
x <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156737/suppl/GSE156737_WT_DMSO.ns.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))


### Get MCF-7
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5557nnn/GSM5557748/suppl/GSM5557748_Si_Ctrl_PRO_Seq_for_METTL3_KD_Rep1_neg_test_noNorm.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5557nnn/GSM5557748/suppl/GSM5557748_Si_Ctrl_PRO_Seq_for_METTL3_KD_Rep1_pos_test_noNorm.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))


# A549
x <- "https://www.encodeproject.org/files/ENCFF275NOU/@@download/ENCFF275NOU.bigWig"
df <-  gsub(x, pattern = ".+download/(ENC.+)$", replacement = "tmp/\\1")
if (! file.exists(df)) download.file(x, destfile = df)
x <- "https://www.encodeproject.org/files/ENCFF979GYA/@@download/ENCFF979GYA.bigWig"
df <-  gsub(x, pattern = ".+download/(ENC.+)$", replacement = "tmp/\\1")
if (! file.exists(df)) download.file(x, destfile = df)


# HFF
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5623nnn/GSM5623305/suppl/GSM5623305_2018-08-07-HFF-hg38-SIC-FW.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5623nnn/GSM5623305/suppl/GSM5623305_2018-08-07-HFF-hg38-SIC-RV.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)


# KELLY
x <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179916/suppl/GSE179916_PROseq_WT_F.bigWig"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)
x <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179916/suppl/GSE179916_PROseq_WT_R.bigWig"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)


# K562
x <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181161/suppl/GSE181161_PROseq_K562_hg38_pl.bigWig"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)
x <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181161/suppl/GSE181161_PROseq_K562_hg38_mn.bigWig"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)


# HEK293
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5401nnn/GSM5401690/suppl/GSM5401690_Proseq.fwd.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5401nnn/GSM5401690/suppl/GSM5401690_Proseq.rev.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)


### DLD1
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4296nnn/GSM4296321/suppl/GSM4296321_PRO-seq-08-Parental-0h.plus.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4296nnn/GSM4296321/suppl/GSM4296321_PRO-seq-08-Parental-0h.minus.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))


# H9 hESCs
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4214nnn/GSM4214080/suppl/GSM4214080_H9_DMSO_rep1_PE1_plus.bigWig"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4214nnn/GSM4214080/suppl/GSM4214080_H9_DMSO_rep1_PE1_minus.bigWig"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
if (! file.exists(df)) download.file(x, destfile = df)


## IMR90
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2828nnn/GSM2828779/suppl/GSM2828779_PD20_PROseq_combined_plus.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))
x <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2828nnn/GSM2828779/suppl/GSM2828779_PD20_PROseq_combined_minus.bw"
df <- gsub(x, pattern = ".+/suppl/", replacement = "") %>% file.path("tmp", .)
df2 <- gsub(df, pattern = "\\.bw", replacement = ".hg38")
if (! file.exists(df)) download.file(x, destfile = df)
system(paste0("~/miniconda3/condabin/conda run -n crossmap CrossMap.py bigwig tmp/hg19_hg38.chain.gz ", df, " ", df2))


## Get pause index ##

library(tidyverse)
library(BRGenomics)
library(rtracklayer)
library(magrittr)
load("tmp/for_pause.rda")

txs <- GenomicFeatures::transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
gb <- genebodies(txs, 300, -300, min.window = 400)
txs <- subset(txs, tx_name %in% gb$tx_name)
pr <- promoters(txs, 0, 100)
gts <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
  AnnotationDbi::select(
    ., keys = AnnotationDbi::keys(.), columns = c("TXNAME", "SYMBOL")
  )

files <- tribble(
  ~cell, ~strand, ~destfile,
  "A549", "+", "tmp/ENCFF979GYA.bigWig",
  "A549", "-", "tmp/ENCFF275NOU.bigWig",
  "HCT116", "+", "tmp/GSE156737_WT_DMSO.ps.hg38.bw.bw",
  "HCT116", "-", "tmp/GSE156737_WT_DMSO.ns.hg38.bw",
  "HeLa", "+", "tmp/GSM5420820_Async_repA_PROseq.pos.bw",
  "HeLa", "-", "tmp/GSM5420820_Async_repA_PROseq.neg.bw",
  "MCF7", "+", "tmp/GSM5557748_Si_Ctrl_PRO_Seq_for_METTL3_KD_Rep1_pos_test_noNorm.hg38.bw",
  "MCF7", "-", "tmp/GSM5557748_Si_Ctrl_PRO_Seq_for_METTL3_KD_Rep1_neg_test_noNorm.hg38.bw",
  "HFF", "+", "tmp/GSM5623305_2018-08-07-HFF-hg38-SIC-FW.bw",
  "HFF", "-", "tmp/GSM5623305_2018-08-07-HFF-hg38-SIC-RV.bw",
  "KELLY", "+", "tmp/GSE179916_PROseq_WT_F.bigWig",
  "KELLY", "-", "tmp/GSE179916_PROseq_WT_R.bigWig",
  "K562", "+", "tmp/GSE181161_PROseq_K562_hg38_pl.bigWig",
  "K562", "-", "tmp/GSE181161_PROseq_K562_hg38_mn.bigWig",
  "HEK293", "+", "tmp/GSM5401690_Proseq.fwd.bw",
  "HEK293", "-", "tmp/GSM5401690_Proseq.rev.bw",
  "B2B", "+", "tmp/GSM5103029_B2BPROseq_V1.pos.bw",
  "B2B", "-", "tmp/GSM5103029_B2BPROseq_V1.neg.bw",
  "DLD1", "+", "tmp/GSM4296321_PRO-seq-08-Parental-0h.plus.hg38.bw",
  "DLD1", "-", "tmp/GSM4296321_PRO-seq-08-Parental-0h.minus.hg38.bw",
  "H9", "+", "tmp/GSM4214080_H9_DMSO_rep1_PE1_plus.bigWig",
  "H9", "-", "tmp/GSM4214080_H9_DMSO_rep1_PE1_minus.bigWig",
  "IMR90", "+", "tmp/GSM2828779_PD20_PROseq_combined_plus.hg38.bw",
  "IMR90", "-", "tmp/GSM2828779_PD20_PROseq_combined_minus.hg38.bw"
) %>% 
  pivot_wider(id_cols = cell, names_from = strand, values_from = destfile)

dir.create("tmp/pidxs", showWarnings = FALSE)

pres <- lapply(
  seq(nrow(files)), function(i) {
    message(i)
    outfile <- paste0("tmp/pidxs/", files$cell[i], ".rds")
    if (! file.exists(outfile)) {
      pos <- import.bw(files$`+`[i])
      strand(pos) <- "+"
      neg <- import.bw(files$`-`[i])
      neg$score <- if(any(neg$score < 0)) -1 * neg$score else neg$score
      strand(neg) <- "-"
      seqs <- c(pos, neg)
      pidx <- getPausingIndices(seqs, pr, gb)
      saveRDS(pidx, outfile)
    } else {
      pidx <- readRDS(outfile)
    }
    
    pidxtbl <- tibble(
      "TXNAME" = gb$tx_name, 
      "pauseIndex" = pidx
    ) %>%
      filter(! is.na(pauseIndex), is.finite(pauseIndex), pauseIndex > 0) %>%
      mutate(TXNAME = gsub(TXNAME, pattern = "\\..+", replacement = "")) %>%
      inner_join(gts, by = "TXNAME") 
    tssrl <- unique(oltxsum$rlregion[oltxsum$type == "TSS"])
    restbltss <- restbl %>%
      filter(rlregion %in% tssrl) %>%
      mutate(cond = case_when(
        log2FoldChange > 0 & padj < 0.05 ~ "S9.6-specific",
        log2FoldChange < 0 & padj < 0.05 ~ "dRNH-specific",
        TRUE ~ "None"
      )) %>%
      filter(cond != "None") %>%
      group_by(cond) %>% 
      slice_min(n = 500, order_by = padj) %>% 
      dplyr::select(cond, SYMBOL=allGenes) %>%
      mutate(SYMBOL = strsplit(SYMBOL, ",")) %>%
      unnest(cols = "SYMBOL")
    inner_join(restbltss, pidxtbl, by = "SYMBOL") %>% 
      mutate(
        cell=files$cell[i]
      )
  }
) %>% bind_rows()


# save(oltxsum, restbl, file = "tmp/for_pause.rda")
ip_cols <- c("dRNH" = "#c386c4", "S9.6"="#82d0e8")

## Main plot of cell lines in both dRNH and S9.6
mainplt <- c("HEK293", "HeLa", "K562")
plt <- pres %>%
  group_by(cell, cond, SYMBOL) %>% 
  summarise(pauseIndex = median(pauseIndex)) %>% 
  filter(cell %in% mainplt) %T>% 
  {
    group_by(., cell, cond) %>% summarise(med=median(pauseIndex)) %>% 
      pivot_wider(id_cols = cell, names_from = cond, values_from = med) %>% 
      mutate(dif = `dRNH-specific`/`S9.6-specific`) %>% 
      arrange(desc(dif)) %>% print()
  } %>%
  ggplot(aes(x = cond, y = pauseIndex, fill = cond)) +
  geom_violin(alpha = .5, width=.9) +
  geom_boxplot(width = .5) +
  ylab("Pause index (log scale)") +
  xlab(NULL) +
  ggpubr::stat_compare_means(
    comparisons = list(c("dRNH-specific", "S9.6-specific")), 
    size = 4.5
  ) +
  scale_y_log10(limits = c(1E-2, 2E5)) +
  scale_fill_manual(values = setNames(ip_cols, nm = paste0(names(ip_cols), "-specific"))) +
  theme_bw(base_size = 16) +
  ggpubr::rremove("legend") +
  facet_wrap(~cell) +
  ggpubr::rotate_x_text(30)
plt
ggsave(plt, filename = "results/Figure_5/pause_3cells_index.png", height = 5, width = 13)

## Supplemental with all other cell lines
plt <- pres %>%
  group_by(cell, cond, SYMBOL) %>% 
  summarise(pauseIndex = median(pauseIndex)) %>% 
  filter(! cell %in% mainplt) %T>% 
  {
    group_by(., cell, cond) %>% summarise(med=median(pauseIndex)) %>% 
      pivot_wider(id_cols = cell, names_from = cond, values_from = med) %>% 
      mutate(dif = `dRNH-specific`/`S9.6-specific`) %>% 
      arrange(desc(dif)) %>% print()
  } %>%
  ggplot(aes(x = cond, y = pauseIndex, fill = cond)) +
  geom_violin(alpha = .5, width=.9) +
  geom_boxplot(width = .5) +
  ylab("Pause index (log scale)") +
  xlab(NULL) +
  ggpubr::stat_compare_means(
    comparisons = list(c("dRNH-specific", "S9.6-specific")), 
    size = 4.5
  ) +
  scale_y_log10(limits = c(1E-2, 1E6)) +
  scale_fill_manual(values = setNames(ip_cols, nm = paste0(names(ip_cols), "-specific"))) +
  theme_bw(base_size = 16) +
  ggpubr::rremove("legend") +
  facet_wrap(~cell) +
  ggtitle("Promoter pausing in genes overlapping dRNH-specific R-loops") +
  ggpubr::rotate_x_text(30)

ggsave(plt, filename = "results/Figure_5/pause_9cells_index.png", height = 10, width = 10)


## Which to use for Genome Browser picture?
topdrnhrl <- restbl %>% 
  slice_min(stat, n = 5)
pres %>% 
  filter(cell == "HEK293", SYMBOL %in% unlist(str_split(topdrnhrl$allGenes, pattern = ","))) %>% 
  arrange(desc(pauseIndex)) %>% 
  slice_max(1)  # ATG3


##### R2R: Size depending on localization
tssrl <- unique(oltxsum$rlregion[oltxsum$type == "TSS"])
locpat <- "(.+):(.+)-(.+):(.+)"
rlgr <- rlregions %>%
  filter(source == "S96") %>% 
  dplyr::select(rlregion, location) %>%
  inner_join(oltxsum) %>% 
  mutate(
    chrom = gsub(location, pattern = locpat, replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = locpat, replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = locpat, replacement = "\\3")),
    width = end - start
  ) %>% dplyr::select(rlregion, width, type)
ggplot(rlgr, aes(x = type, y = width, fill = type)) +
  geom_jitter(width = .15, alpha=.2) +
  scale_y_log10()




##### R2R: False negatives

rlsamples %>% 
  filter(mode == "MapR", genome == "mm10") %>%
  mutate(
    condition = case_when(
      grepl(condition, pattern = "\\-t[0-9]{1}$") ~ "RNase T",
      condition == "RNaseA" ~ "RNase A",
      grepl(other, pattern = "ActD") ~ "Actinomycin D",
      condition == "Input" & mode == "MapR" ~ "MNase", # This one was mislabeled
      condition == "pAMNase" & mode == "MapR" ~ "MNase",
      condition == "RNH" ~ "RNase H1",
      condition == "RNHdMNase" ~ "NBSS",  # This one was mislabeled
      TRUE ~ condition
    )
  ) %>% 
  group_by(condition) %>% 
  mutate(url = file.path(RLSeq:::RLBASE_URL, coverage_s3)) %>%
  select(rlsample, label, condition, prediction, mode, tissue, other, url) %>% View()

# Wrangle and repair condition names for plotting
conds <- rlsamples %>% 
  filter(label == "NEG", prediction == "POS") %>% 
  mutate(
    condition = case_when(
      grepl(condition, pattern = "\\-t[0-9]{1}$") ~ "RNase T",
      condition == "RNaseA" ~ "RNase A",
      grepl(other, pattern = "ActD") ~ "Actinomycin D",
      condition == "Input" & mode == "MapR" ~ "MNase", # This one was mislabeled
      condition == "pAMNase" & mode == "MapR" ~ "MNase",
      condition == "RNH" ~ "RNase H1",
      condition == "RNHdMNase" ~ "NBSS",  # This one was mislabeled
      TRUE ~ condition
    )
  ) %>% 
  group_by(condition) %>% 
  mutate(url = file.path(RLSeq:::RLBASE_URL, coverage_s3)) %>%
  select(rlsample, label, condition, mode, tissue, other, url) 

## Meta RLFS for Input
toplt <- conds %>%
  filter(condition == "Input") %>%
  pull(rlsample) %>% 
  lapply(function(x) {
    object <- RLSeq::RLRangesFromRLBase(x)
    rlfsRes <- rlresult(object, resultName = "rlfsRes")
    tibble(
      shift = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifts,
      zs = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores,
      rlsample = x
    )
  })

toplt %>%
  bind_rows() %>% 
  group_by(rlsample) %>% 
  mutate(zs = minmax(zs)) %>% 
  ggplot(aes(x = shift, y = zs, color = rlsample, group=rlsample)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  geom_line(size=1.5, alpha = .8) +
  ylab("Z score (scaled)") +
  xlab("Distance from RLFS (bp)") +
  theme_prism(border = TRUE, base_size = 16) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  ggtitle("Input control RLFS analysis (DRIP)")


## Meta RLFS for MNase
toplt <- conds %>%
  filter(condition == "MNase") %>%
  pull(rlsample) %>% 
  lapply(function(x) {
    object <- RLSeq::RLRangesFromRLBase(x)
    rlfsRes <- rlresult(object, resultName = "rlfsRes")
    tibble(
      shift = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifts,
      zs = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores,
      rlsample = x
    )
  })

toplt %>%
  bind_rows() %>% 
  group_by(rlsample) %>% 
  mutate(zs = minmax(zs)) %>% 
  ggplot(aes(x = shift, y = zs, color = rlsample, group=rlsample)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  geom_line(size=1.5, alpha = .8) +
  ylab("Z score (scaled)") +
  xlab("Distance from RLFS (bp)") +
  theme_prism(border = TRUE, base_size = 16) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  ggtitle("MNase control RLFS analysis (MapR)")


## Put mouseRLFS on genome browser
infile <- "https://rlbase-data.s3.amazonaws.com/rlfs-beds/mm10.rlfs.bed"
outfile <- "tmp/mm10.rlfs.bb"
if (! file.exists(outfile)) {
  valr::read_bed12(infile) %>% 
    mutate(score = 1, strand = ".", name = paste0("RLFS_", row_number())) %>% 
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = T, 
      seqinfo =  GenomeInfoDb::getChromInfoFromUCSC("mm10", as.Seqinfo = TRUE)
    ) %>% 
    rtracklayer::export.bb(con = outfile)
  aws.s3::put_object(file = outfile, object = paste0("misc/mm10.rlfs.bb"), bucket = RLSeq:::RLBASE_S3, multipart = TRUE, show_progress = TRUE)
}

conds %>%
  filter(condition == "MNase") %>% 
  pull(url)


## Donut chart to summarize false neg
plt <- plot_ly(type = "pie", size = I(.01)) %>%
  add_pie(
    data=tally(conds),
    values = ~n, labels = ~condition, textinfo='label+value',
    insidetextorientation='horizontal', hole=.6, rotation=186
  ) %>%
  layout(showlegend = FALSE, margin = list(l = 100, r = 100, t=100, b=100),
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
plt
save_image(plt, file = "results/Figure_3/false_neg_cond_donut.svg", width = 300, height = 300, format = "svg")




##### R2R: reproducibility with down sample

## Get differential RL Regions
restbls <- pbapply::pblapply(seq(10), function(i) {
  
  # Make random
  set.seed(i)
  
  # Get the colData
  cd <- ctsPos@colData %>% 
    as_tibble() %>%
    inner_join(
      select(rlsamples, rlsample, ip_type),
      by = c("experiment" = "rlsample")
    ) 
  
  # Number of studies to randomly select from S9.6 data
  num_to_sel <- rlsamples %>% 
    filter(rlsample %in% cd$experiment & ip_type == "dRNH") %>% 
    pull(study) %>% unique() %>% length()
  
  # Studies to select from
  studies_selected <- rlsamples %>% 
    filter(rlsample %in% cd$experiment & ip_type == "S9.6") %>% 
    pull(study) %>%
    unique() %>%
    sample(size = num_to_sel)
  
  # Get the S9.6 sample IDs to include (also down sample to hit N=56)
  sampl <- rlsamples %>% 
    filter(
      rlsample %in% cd$experiment,
      ip_type == "S9.6",
      study %in% studies_selected
    ) %>% 
    slice_sample(n=56) %>% 
    pull(rlsample)    
  
  # Wrangle DDS and run analysis
  drns <- cd$experiment[cd$ip_type == "dRNH"]
  ctsPossmall <- ctsPos[,ctsPos$experiment %in% c(sampl, drns)]
  matsmall <- ctsPossmall@assays@data$cts
  cdsmall <- cd[cd$experiment %in%  c(sampl, drns),]
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = matsmall[,cdsmall$experiment], colData = cdsmall, design = ~ip_type
  )
  dds <- DESeq2::DESeq(dds)
  DESeq2::results(dds, contrast = c("ip_type", "S9.6", "dRNH")) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "rlregion") %>%
    as_tibble() %>%
    mutate(rlregion = gsub(rlregion, pattern = "All_", replacement = "")) %>% 
    inner_join(rlregions, by = "rlregion") %>%
    arrange(padj) %>% 
    mutate(
      section = paste0("permutation_", i),
      n_section = length(sampl)
    )
})


df <- restbl %>% 
  mutate(
    section = "Full",
    n_section = 191
  ) %>% 
  bind_rows(bind_rows(restbls))
annodf <- df %>% 
  select(section, n_selected=n_section) %>%
  distinct() %>% 
  column_to_rownames("section")
corr <- df %>% 
  select(rlregion, stat, section) %>% 
  pivot_wider(id_cols = rlregion, names_from = section, values_from = stat) %>% 
  column_to_rownames("rlregion") %>% 
  cor() 

# Sets the minimum (0), the maximum (15), and the increasing steps (+1) for the color scale
# Note: if some of your genes are outside of this range, they will appear white on the heatmap
library(pheatmap)
library(RColorBrewer)

breaksList = seq(0, 1, by = .001)

# Plots the first heatmap
pheatmap(
  corr, 
  color = colorRampPalette((brewer.pal(n = 7, name = "BrBG")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
  breaks = breaksList,
  cluster_cols = F,
  cluster_rows = F,
  angle_col = 45,
  fontsize = 14,
  display_numbers = T, number_color = "white", fontsize_number = 18,
  main = "Correlation of differential analysis between permutations\nPearson correlation of Wald statistic"
) 

# 4-way plot
df %>% 
  select(rlregion, stat, section) %>% 
  pivot_wider(id_cols = rlregion, names_from = section, values_from = stat) %>% 
  ggplot(aes(x = Full, y = permutation_1)) +
  geom_hline(yintercept = 0, color="red", linetype="dashed") +
  geom_vline(xintercept = 0, color="red", linetype="dashed") +
  geom_point() +
  theme_gray(16) +
  stat_smooth(method = "lm", formula = y ~ 0 + x) +
  ggtitle(
    "Downsampled vs full dataset",
    subtitle = "Differential analysis results (comparison of Wald statistic)"
  ) +
  xlab("Full dataset (S9.6 vs dRNH Wald stat)") +
  ylab("Downsampled dataset (S9.6 vs dRNH Wald stat)") +
  ggpmisc::stat_poly_eq(
    formula = y ~ 0 + x,
    size=5,
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    parse = TRUE,
  )



###### Percent of genome covered by consensus regions
locpat <- "(.+):(.+)-(.+):(.+)"

### Total coverage
rlgr <- rlregions %>%
  dplyr::select(location) %>%
  mutate(
    chrom = gsub(location, pattern = locpat, replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = locpat, replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = locpat, replacement = "\\3"))
  ) %>%
  dplyr::select(-location) %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%
  GenomicRanges::reduce() 

rs <- rlgr %>% 
  GenomicRanges::width() %>%
  sum()

wg <- 3209286105
prop <- rs/wg  # From http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics
prop


### Shared coverage
rlgr <- rlregions %>%
  filter(source == "dRNH S96") %>% 
  dplyr::select(location) %>%
  mutate(
    chrom = gsub(location, pattern = locpat, replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = locpat, replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = locpat, replacement = "\\3"))
  ) %>%
  dplyr::select(-location) %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%
  GenomicRanges::reduce() 

rs <- rlgr %>% 
  GenomicRanges::width() %>%
  sum()

wg <- 3209286105
prop <- rs/wg  # From http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics
prop


### PCT plot
plt <- plot_ly(type = "pie") %>%
  add_pie(data=tibble(
    group=c("  ", " "),
    size=c(rs, (wg - rs))
  ) %>% mutate(pct = round(100*size/sum(size), 2)),
  values = ~pct, labels = ~group, textinfo='label+value',
  insidetextorientation='horizontal', hole=.6, rotation=186) %>%
  layout(showlegend = FALSE, margin = list(l = 100, r = 100, t=100, b=100),
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
save_image(plt, file = "results/Rloop_genome_donut.svg", format = "svg")

lapply(cons2, function(x) {
  x %>% 
    mutate(width = end - start) %>% 
    pull(width) %>%
    {as_tibble(sum(.)/1E6)}
}) %>% 
  bind_rows(.id = "sel") %>% 
  filter(sel != "S9.6") %>% 
  mutate(
    group = gsub(sel, pattern = "\\-.+", replacement = "")
  ) %>% 
  group_by(group) %>% 
  {setNames(group_split(.), nm = group_keys(.)[[1]])} %>% 
  lapply(
    function(x) {
      tibble(
        sel = "Other RL-Regions", value = (rs/1E6) - sum(x$value), group = x$group[1]
      ) %>% bind_rows(x)
    } 
  ) %>% 
  bind_rows() %>% 
  mutate(sel = factor(sel),
         sel= relevel(sel, ref="Other RL-Regions"),
         group = factor(group, levels = rev(c("dRNH", "S9.6")))) %>% 
  ggplot(aes(x = group, y = value, fill=sel)) +
  geom_col(color="black") +
  coord_flip() +
  xlab(NULL) +
  ylab("R-loop region genomic coverage (Mb)") +
  theme_bw(base_size = 20, base_line_size = .5) +
  scale_fill_manual(
    values = c(
      "dRNH-only" = "#d7a5cb",
      "dRNH-shared" = "#997691",
      "S9.6-only" = "#83cfe7",
      "S9.6-shared" = "#81aeb8",
      "Other RL-Regions" = "#d6d6d6"
    )
  ) +
  guides(fill=guide_legend(title = NULL)) -> plt
plt
ggsave(
  plt, filename = "results/Figure_4/prop_RLRegion_occupied.png",
  width=851, height=511, units = "px", dpi = 96
)






