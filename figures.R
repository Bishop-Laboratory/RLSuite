### Main script for generating figures (eventually should be notebook) ###
library(reticulate)
library(tidyverse)
library(kableExtra)
library(tidymodels)
library(GetoptLong)
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

## Table 1 -- univariate analysis ##
rlsamples %>%
  tableone::CreateTableOne(
    data = .,
    test = FALSE, 
    vars = c("genome", "mode", "ip_type", "label")
  ) %>%
  tableone::kableone() 

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


# Num peaks across conditions 
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
  filter(! is.na(numPeaks)) %>%
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

### Figure 2 ###
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

# With just selected samples
pos_pos <- c(
  "SRX3084731",
  "SRX4776642",
  "SRX1070677",
  "SRX7583981",
  "SRX2683608",
  "SRX7347914"
)
pos_neg <- c(
  "SRX10229652",
  "SRX5696405",
  "SRX2642969",
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


## Metaplots of RLFS (main)
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
ggsave(plt, filename = "results/Figure_1/rlfs_zscore_plot.svg", height = 5, width = 7.5)

rlsamples %>% group_by(label, prediction) %>% tally()

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
chm <- ComplexHeatmap::pheatmap(heatdata, color = colors1, main = "Percentage samples labels by Mode",
                                fontsize = 13, cluster_rows = FALSE, cluster_cols = FALSE,
                                display_numbers = TRUE, fontsize_number = 10.5,
                                annotation_col = annocol, 
                                show_colnames = FALSE, name = "Representation (%)",
                                annotation_colors = cat_cols, use_raster=FALSE,
                                annotation_row = annorow[rownames(heatdata),,drop=FALSE],
                                number_color = "black", angle_col = "0")

pdf(qq("results/Figure_1/heatmap_representation.pdf"), 
    width = 8, height = 8)
ComplexHeatmap::draw(chm)
dev.off()



### Figure 3/S3
dir.create("results/Figure_3", showWarnings = FALSE)

## Genome Browser -- get scaling factors for bigWigs
## From https://www.biostars.org/p/413626/#414440
rlcts <- RLHub::rlregions_counts()
cts_mat <- rlcts@assays@data$cts
normFact <- edgeR::calcNormFactors(object = cts_mat, method="TMM")
libSize <- colSums(cts_mat)
sizeFactors <- normFact * libSize / 1000000
sizeFactors.Reciprocal <- 1/sizeFactors
sizeFactTbl <- tibble(
  "rlsample" = names(sizeFactors),
  "sizeFactors" = sizeFactors,
  "sizeFactors.Reciprocal" = sizeFactors.Reciprocal
) %>% write_csv(file = "data/sizeFactors.csv")

# Re-run deepTools with the scaling (done from the CLI)
# 1. Predict fragment length
# macs3 predictd -i "../RLBase-data/rlbase-data/rlpipes-out/bam/SRX2455193/SRX2455193_hg38.bam" -g hs
tribble(
  ~rlsample, ~predict_fragsize,
  "ERX3974964", 245,
  "ERX3974965", 255,
  "SRX2455189", 197,
  "SRX2455193", 182,
  "SRX1025894", 191,
  "SRX1025896", 294,
  "SRX2683605", 296,
  "SRX2675009", 260
)
# 2. Run deeptools with scaling factors and extended fragments
# bamCoverage -b "bam/ERX3974964/ERX3974964_hg38.bam" -o coverage_scaled/ERX3974964_hg38.scale.bw --scaleFactor 0.07312689 -p 44 -e 245 --ignoreDuplicates
# bamCoverage -b "bam/ERX3974965/ERX3974965_hg38.bam" -o coverage_scaled/ERX3974965 --scaleFactor 0.2667499 -p 44 -e 255 --ignoreDuplicates
# bamCoverage -b "bam/SRX2455189/SRX2455189_hg38.bam" -o coverage_scaled/SRX2455189_hg38.scaled.bw --scaleFactor 0.7608039 -p 44 -e 197 --ignoreDuplicates
# bamCoverage -b "bam/SRX2455193/SRX2455193_hg38.bam" -o coverage_scaled/SRX2455193_hg38.scaled.bw --scaleFactor 0.5406119 -p 44 -e 182 --ignoreDuplicates
# bamCoverage -b "bam/SRX1025894/SRX1025894_hg38.bam" -o coverage_scaled/SRX1025894_hg38.scale.bw --scaleFactor 0.1704441 -p 44 -e 191 --ignoreDuplicates
# bamCoverage -b "bam/SRX1025896/SRX1025896_hg38.bam" -o coverage_scaled/SRX1025896_hg38.scale.bw --scaleFactor 0.7183404 -p 44 -e 294 --ignoreDuplicates
# bamCoverage -b "bam/SRX2683605/SRX2683605_hg38.bam" -o coverage_scaled/SRX2683605_hg38.scale.bw --scaleFactor 0.6507813 -p 44 -e 296 --ignoreDuplicates
# bamCoverage -b "bam/SRX2675009/SRX2675009_hg38.bam" -o coverage_scaled/SRX2675009_hg38.scale.bw --scaleFactor 0.9435923 -p 44 -e 260 --ignoreDuplicates


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
pltdat <- RLSeq::plotEnrichment(rlr, returnData = TRUE, 
                                pred_POS_only = FALSE,
                                rlbaseRes = sampleAnno, 
                                rlsamples = rlsamples,
                                label_POS_only = FALSE) 

# Genes
plt <- plot_multi_feature(features = unique(pltdat$Transcript_Features$type), db = "Transcript_Features", 
                          factororder = c("TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS"),
                          pltdat = pltdat,
                          lmts = c(-6, 12),
                          axislmts = c(-6, 15),
                          yb = 11.5) +
  ggtitle("Genic Feature Enrichment in RLBase samples", subtitle = NULL) + xlab(NULL)
ggsave(plt, filename = "results/Figure_3/genic_enrichment.svg", height = 5, width = 11)

# G4Q
pltdat2 <- pltdat
pltdat2$G4Qpred$type <- "G4Q Predicted"
pltdat2$G4Qexp$type <- "G4Q ChIP"
feats <- c("G4Q Predicted", "G4Q ChIP")
plt <- plot_multi_feature(features = feats, db = c("G4Qpred", "G4Qexp"), 
                          factororder = feats,
                          pltdat = pltdat2,
                          lmts = c(-6, 7),
                          axislmts = c(-6, 8),
                          yb = 11.5) +
  ggtitle("G4Q Enrichment in RLBase samples", subtitle = NULL) + xlab(NULL)
ggsave(plt, filename = "results/Figure_3/g4q_enrichment.svg", height = 5, width = 7)


### Figure 4/S4 ###
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
) %>% GenomicRanges::reduce()

# Density plot of peak sizes
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
# computeMatrix scale-regions -S ../RLBase-data/rlbase-data/rlregions/rlregions_S96.bw ../RLBase-data/rlbase-data/rlregions/rlregions_dRNH.bw -R ~/.rlpipes_genomes/hg38/hg38.ensGene.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 -o matrix.mat.gz -p 44 --verbose
metap <- read_tsv("matrix.mat.gz", skip = 1, col_names = FALSE)
metadat <- read_lines("matrix.mat.gz", n_max = 1)
metadat <- gsub(metadat, pattern = "@", replacement = "") %>% jsonlite::fromJSON()
mat <- metap %>%
  select(-c(1:6)) 
mat %>%
  colMeans(na.rm = TRUE) -> cm
divs <- setNames(c(42, 184), nm = c("rlregions_dRNH", "rlregions_S96"))
pltdats <- tibble("score" = cm,
              "loc" = rep(seq(1100), 2),
              "div" = c(rep(divs[metadat$sample_labels[1]], metadat$sample_boundaries[2]),
                        rep(divs[metadat$sample_labels[2]], metadat$sample_boundaries[2])),
              "label" = c(rep(metadat$sample_labels[1], metadat$sample_boundaries[2]),
                          rep(metadat$sample_labels[2], metadat$sample_boundaries[2]))) %>%
  # group_by(label) %>%
  mutate(
    # score = scales::rescale(x = score, to = c(0, 1)),
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
      pct_cons = 100 * ((score / 10) / 42),  # only 42 samples in the dRNH group ,
      source = "dRNH"
    )
  } else {
    x <- mutate(
      x, 
      pct_cons = 100 * ((score / 10) / 184),  # Compared to 184 for S9.6
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
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  coord_flip() +
  ylab("Peak location (%)") +
  xlab(NULL) +
  scale_fill_manual(values = txfeat_cols, guide = guide_legend(title = "Feature", reverse = TRUE))
ggsave(plt, filename = "results/Figure_4/consensus__txplot_2.svg")

# Check RLFS
ol <- ChIPpeakAnno::findOverlapsOfPeaks(
  GenomicRanges::makeGRangesFromDataFrame(cons$dRNH),
  GenomicRanges::makeGRangesFromDataFrame(cons$S9.6),
  RLSeq:::getRLFSAnno(object = rlr)
)
pdf("results/Figure_4/RLFS_venn.pdf", height = 4, width = 7)
ChIPpeakAnno::makeVennDiagram(ol, fill = c("#e2a3e3", "#82d0e8", "white"), 
                              NameOfPeaks = c("dRNH", "S9.6", "RLFS"))
dev.off()


### Figure 5/S5 ###

## PCA from RLRegion counts
ip_correct <- rlsamples$rlsample[rlsamples$ip_type %in% c("S9.6", "dRNH")]
ctsPos <- cts[
  rownames(cts) %in% rlregions$rlregion[rlregions$source == "dRNH S96"] 
  ,
  ! is.na(cts$prediction) &
    cts$prediction == "POS" & 
    cts$experiment %in% ip_correct &
    cts$label == "POS" &
    cts$numPeaks > 5000]
ctsmat <- ctsPos@assays@data$cts
vstmat <- DESeq2::vst(ctsmat)
pcdata <- prcomp(t(vstmat))
pcsmall <- pcdata$x[,1:7]
# tsdat <- Rtsne::Rtsne(pcsmall, dims = 2, pca = FALSE, theta = 0)
# tdat <- as.data.frame(tsdat$Y)
# rownames(tdat) <- rownames(pcdata$x)  
# colnames(tdat) <- c("TSNE_1", "TSNE_2")
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
  geom_point() +
  # scale_color_manual(values = setNames(RLSeq:::auxdata$mode_cols$col,
  #                                      nm = RLSeq:::auxdata$mode_cols$mode)) +
  scale_color_manual(values = ip_cols) +
  theme_bw(base_size = 14) 
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
restbl <- DESeq2::results(dds) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "rlregion") %>%
  as_tibble() %>%
  inner_join(rlregions, by = "rlregion") %>%
  arrange(padj)

## Volcano plot
volcplt <- restbl %>%
  EnhancedVolcano::EnhancedVolcano(
    lab = gsub(.$rlregion, pattern = "All_", replacement = ""), x = "log2FoldChange", y = "padj", 
    pCutoff = 1E-50, FCcutoff = 2,labSize = 4,
    title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none"
  ) +
  theme_bw(base_size = 14) +
  ggpubr::rremove("legend")
ggsave(plot = volcplt, filename = "results/Figure_5/volcano.svg", height = 5, width = 6)
volcplt +
  theme_void() +
  ggpubr::rremove("legend") -> volcplt2
ggsave(plot = volcplt2, filename = "results/Figure_5/volcano.png", height = 5, width = 6)

# Gene set enrichment
library(enrichR)
ol <- restbl %>%
  mutate(group = ifelse(log2FoldChange > 0, "Over-abundant (S9.6 vs dRNH)", 
                        ifelse(log2FoldChange < 0, "Under-abundant (S9.6 vs dRNH)", "Nothing"))) %>%
  filter(padj < .05) %>%
  group_by(group) %>%
  slice_min(order_by = padj, n = 1000) %>%
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
  slice_max(order_by = combined_score, n = 7) %>%
  pull(Term)
barplt <- pltdata %>%
  filter(Term %in% terms) %>%
  mutate(combined_score = ifelse(! is.na(combined_score), ifelse(combined_score < 1, 0, combined_score), 0)) %>%
  arrange(ip_type, (combined_score)) %>%
  mutate(
    Term = gsub(Term, pattern = " \\(GO.+", replacement = ""),
    Term = stringr::str_wrap(Term, width = 45),
    Term = factor(Term, levels = unique(Term)),
    ip_type = gsub(ip_type, pattern = "S96", replacement = "S9.6")
  ) %>%
  ggplot(aes(x = Term, y = combined_score, fill = ip_type)) +
  geom_col(position = position_dodge(.8)) +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab("Combined Score") +
  xlab(NULL) +
  scale_fill_manual(values = ip_cols)
ggsave(barplt, filename = "results/Figure_5/barplt.svg", height = 8, width = 8)

# Feature plots
to_plt <- restbl$rlregion[1:2]

plts <- lapply(to_plt, function(pt) {
  plt <- tibble(
    "cts" = vstmat[pt,],
    "name" = names(vstmat[pt,] )
  ) %>% inner_join(pltdat) %>%
    ggplot(aes(x = PC1, y = PC2, color = cts, group = tissue)) +
    geom_hline(yintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
    geom_vline(xintercept = 0, color="grey", alpha = 1, linetype = "dashed") +
    geom_point() +
    scale_color_viridis_b(option = "A", direction = -1) +
    theme_bw(base_size = 14) 
  ggsave(plot = plt, filename = paste0("results/Figure_5/pca_feature_plot_", pt,".svg"), height = 5, width = 7)
})

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

# Show this on genome browser
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


### Figure 6: Conserved vs Variable ###

dir.create("results/Figure_6", showWarnings = FALSE)

## Define conserved vs variable
# Due to a mistake in the conservation score, we need to use a scaling factor
# and convert it. 61 divisor for dRNH and 270 divisor for S9.6
# Original (flawed) calc divided number (N) of samples overlapping with 
# windows by the divisor instead of the proper scaling factor (42, 184)
# To reverse this, simply multiply by the divisor / proper factor
scale_fact <- c("S9.6" = 270/184, "dRNH" = 61/42)
# Apply scaling factor
rlregions_cons <- rlregions %>%
  mutate(
    pct_cons = case_when(
      source == "dRNH" ~ conservation_score * scale_fact['dRNH'],
      source == "S9.6" ~ conservation_score * scale_fact['S9.6'],
      TRUE ~ conservation_score * mean(scale_fact)
    )
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
ggsave(plt, filename = "results/Figure_6/pct_cons_ranks.svg", height = 6, width = 9)

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
ggsave(plt, filename = "results/Figure_6/pct_cons_ridges.svg")

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
                    guide = guide_legend(title = "Feature", reverse = TRUE)) +
  scale_x_discrete(limits = rev)
ggsave(plt, filename = "results/Figure_6/pct_cons_txfeats.svg")

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
  filter(allGenes %in% xrn2genes) %>%
  View()

num_sel <- 5
db <- "KEGG_2021_Human"
# db <- "MSigDB_Hallmark_2020"
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
  mutate(combined_score = ifelse(! is.na(combined_score),
                                 ifelse(combined_score < 1, 0, 
                                        combined_score), 0)) %>%
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
ggsave(plt, filename = "results/Figure_6/path_enrich_tf.svg", height = 8, width = 6)

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
  filename = "results/Figure_6/housekeeping_genes.svg",
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

## Table S3 Differentially-abundant RL Regions
dir.create("results/Table_S3/", showWarnings = FALSE)
restbl %>%
  mutate(rlregion = gsub(rlregion, pattern = "All_", replacement = "")) %>%
  select(rlregion, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  write_csv("results/Table_S3/differentially_abundant_rlregions.csv")


### Pausing index
library(BRGenomics)
library(rtracklayer)
files <- tribble(
  ~cell, ~strand, ~URL,
  "A549", "+", "https://www.encodeproject.org/files/ENCFF979GYA/@@download/ENCFF979GYA.bigWig",
  "A549", "-", "https://www.encodeproject.org/files/ENCFF275NOU/@@download/ENCFF275NOU.bigWig"
) %>% mutate(destfile = gsub(URL, pattern = ".+download/(ENC.+)$", replacement = "tmp/\\1"))
dir.create("tmp", showWarnings = FALSE)
sapply(files$URL, function(x) {
  destfile <- gsub(x, pattern = ".+download/(ENC.+)$", replacement = "tmp/\\1")  
  download.file(x, destfile = destfile)
})
pos <- import.bw(files$destfile[1])
strand(pos) <- "+"
neg <- import.bw(files$destfile[2])
neg$score <- -1 * neg$score
strand(neg) <- "-"
seqs <- c(pos, neg)
txs <- GenomicFeatures::transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
gb <- genebodies(txs, 300, -300, min.window = 400)
txs <- subset(txs, tx_name %in% gb$tx_name)
pr <- promoters(txs, 0, 100)
pidx <- getPausingIndices(seqs, pr, gb)
gts <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 %>%
  AnnotationDbi::select(
    ., keys = AnnotationDbi::keys(.), columns = c("TXNAME", "SYMBOL")
  )
pidxtbl <- tibble(
  "TXNAME" = gb$tx_name, 
  "pauseIndex" = pidx
) %>%
  filter(! is.na(pauseIndex), is.finite(pauseIndex), pauseIndex > 0) %>%
  mutate(TXNAME = gsub(TXNAME, pattern = "\\..+", replacement = "")) %>%
  inner_join(gts, by = "TXNAME") 
tssrl <- unique(oltxsum$rlregion[oltxsum$type == "TSS"])
respid <- restbl %>%
  filter(rlregion %in% tssrl) %>%
  mutate(cond = case_when(
    log2FoldChange > 0 & padj < .05 ~ "S9.6-specific",
    log2FoldChange < 0 & padj < .05 ~ "dRNH-specific",
    TRUE ~ "None"
  )) %>%
  filter(cond != "None") %>%
  dplyr::select(cond, SYMBOL=allGenes) %>%
  mutate(SYMBOL = strsplit(SYMBOL, ",")) %>%
  unnest(cols = "SYMBOL") %>%
  inner_join(pidxtbl, by = "SYMBOL")

library(magrittr)
respid %T>%
  {
    group_by(., cond) %>% summarise(median(pauseIndex)) %>% print()
  } %>%
  ggplot(aes(x = cond, y = pauseIndex, fill = cond)) +
  geom_violin(alpha = .5) +
  geom_boxplot(width = .5) +
  scale_y_log10() +
  ylab("Pause index (log scale)") +
  xlab(NULL) +
  ggpubr::stat_compare_means(comparisons = list(c("dRNH-specific", "S9.6-specific")), 
                             label = "p.signif", size = 6) +
  scale_fill_manual(values = setNames(ip_cols, nm = paste0(names(ip_cols), "-specific"))) +
  theme_bw(base_size = 14) +
  ggpubr::rremove("legend") -> plt
ggsave(plt, filename = "results/Figure_5/pause_index.svg", height = 6, width = 6)


## Percent of genome
locpat
rlgr <- rlregions %>%
  dplyr::select(location) %>%
  mutate(
    chrom = gsub(location, pattern = locpat, replacement = "\\1"),
    start = as.numeric(gsub(location, pattern = locpat, replacement = "\\2")),
    end = as.numeric(gsub(location, pattern = locpat, replacement = "\\3"))
  ) %>%
  dplyr::select(-location) %>%
  GenomicRanges::makeGRangesFromDataFrame() %>% reduce() 

rs <- rlgr %>% 
  width() %>%
  sum()

rs/3209286105  # From http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics


