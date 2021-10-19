### Main script for generating figures (eventually should be notebook) ###
library(reticulate)
library(tidyverse)
library(kableExtra)
library(tidymodels)
library(GetoptLong)
library(VennDiagram)
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

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create(file.path("results", "tables"), showWarnings = FALSE)

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

## RLFS plots
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

# Show accuracy of the model

# Show some of the RLFS analysis results
# Make a metaplot with confidence intervals 
# Split by predictions
rlfsRes <- RLHub::rlfs_res()
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

plt <- toplt %>%
  rename(Label = label, Prediction=prediction) %>%
  ggplot(aes(x = shifts, y = Zscore, color = Label, fill=Label)) +
  geom_vline(xintercept = 0, linetype='dashed', color='grey', alpha=.5) +
  stat_smooth(method="loess", span=.01, se=TRUE, alpha=.3) +
  facet_grid(Prediction ~ Label, labeller = label_both) +
  ylab("Z-score (scaled)") +
  xlab("Distance from RLFS (bp)") +
  theme_prism(border = TRUE, base_size = 10) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  scale_color_manual(values = c("NEG" = "#9e679e", "POS" = "#67889e")) +
  scale_fill_manual(values = c("NEG" = "#9e679e", "POS" = "#67889e")) +
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

# Num peaks across conditions 
plt <- rlsamples %>%
  filter(prediction == "POS" & label == "POS") %>%
  mutate(numPeaks = numPeaks / 1000) %>%
  group_by(mode) %>%
  mutate(numpeakmed = median(numPeaks)) %>%
  ungroup() %>%
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


## Genome Browser -- get scaling factors for bigWigs
## From https://www.biostars.org/p/413626/#414440
rlcts <- RLHub::rlregions_counts()
cts <- rlcts@assays@data$cts
normFact <- edgeR::calcNormFactors(object = cts, method="TMM")
libSize <- colSums(cts)
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

### Figure 2/S2 - Prediction validation ###

dir.create("results/Figure_2", showWarnings = FALSE)

# Num of peaks across condition + RNH status
DODGE_POS <- .8
df_p_val <- filter(rlsamples, ! is.na(prediction)) %>% 
  mutate(numPeaks = numPeaks / 1000) %>%
  rstatix::group_by(mode) %>%
  filter(! is.na(prediction), length(which(prediction == "POS")) > 2,
         length(which(prediction == "NEG")) > 2) %>%
  rstatix::wilcox_test(numPeaks ~ prediction, p.adjust.method = "fdr") %>%
  rstatix::add_significance(p.col = "p") %>% 
  rstatix::add_xy_position(dodge = DODGE_POS, x = "mode")
plt <- rlsamples %>%
  filter(! is.na(prediction)) %>%
  mutate(numPeaks = numPeaks / 1000) %>%
  right_join(df_p_val, by = "mode") %>%
  ggplot(mapping = aes(x = mode,
                       y = numPeaks, fill = prediction)) +
  geom_boxplot(color = "black", width = .7, outlier.shape = NA,
               position = position_dodge(DODGE_POS)) +
  geom_jitter(position = position_jitterdodge(.2), alpha = .5) +
  ylab("Peaks Called (thousands)") +
  xlab(NULL) +
  labs(title = "Peaks called per sample", subtitle = "Split by Mode and Prediction") +
  guides(y = "prism_offset_minor") + 
  theme_prism(base_size = 14) + 
  scale_fill_manual(values = setNames(RLSeq:::auxdata$prediction_cols$col, nm = RLSeq:::auxdata$prediction_cols$prediction)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1)) +
  add_pvalue(df_p_val,  
             xmin = "xmin", 
             xmax = "xmax",
             inherit.aes = F) 
ggsave(plt, filename = "results/Figure_2/peaks_called_prediction.svg")


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
pdf(qq("results/Figure_2/corr_heatmap.pdf"), 
    width = 11, height = 8)
ComplexHeatmap::draw(hm)
dev.off()

## Enrichment in features known for R-loop formation ##

## Build the table of annotation types with the 
## Significance of enrichment (wilcox) within R-loop mapping samples
rlfeatres <- sampleAnno %>%
  left_join(rlsamples, by = c("experiment" = "rlsample")) %>%
  filter(! is.infinite(stat_fisher_rl), genome == "hg38",
         ! is.na(prediction), ! is.na(stat_fisher_rl)) %>%
  group_by(db, type) %>%
  mutate(stat_fisher_rl = (stat_fisher_rl+1) / (stat_fisher_shuf+1)) %>%
  select(c(prediction, stat_fisher_rl)) %>%
  na.omit() %>%
  nest(data = c(prediction, stat_fisher_rl)) %>%
  mutate(
    fit = map(data, ~ wilcox.test(stat_fisher_rl ~ prediction, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(p.value, method = "bonferroni")
  ) %>%
  select(-data, -fit, -method, -alternative) %>%
  arrange(padj)
write_csv(rlfeatres, file.path("results", "tables", "RLBase_anno_enrichment_summarized.csv"))

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
ggsave(plt, filename = "results/Figure_2/genic_enrichment.svg", height = 5, width = 11)

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
ggsave(plt, filename = "results/Figure_2/g4q_enrichment.svg", height = 5, width = 7)



## Compare with RLBPs
rlbpGeneList <- rlbps$geneName[rlbps$geneName %in% humanGenes$SYMBOL & rlbps$combinedScore > 0]
rlfeatres %>%
  filter(
    type %in% humanGenes$SYMBOL
  ) %>%
  mutate(
    RLBP = type %in% {{ rlbpGeneList }},
    padj = -log10(padj),
    statisticrank = rank(-statistic),
    label = case_when(
      RLBP & padj > 30 ~ type, TRUE ~ ""
    )
  ) %>%
  arrange(RLBP) %>%
  ggplot(aes(x = statistic, y = padj, color = RLBP, label=label)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 10, force = 5,
                           nudge_y = -1, max.time = 5, max.iter = 20000, nudge_x = 100) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "#dedede"))  +
  geom_rug(sides = "b") +
  xlab("Wilcox Statistic") + ylab("P Adjusted Value (-log10)") +
  ggtitle("R-loop binding factors vs known RLBPs") +
  theme_bw(base_size = 14) +
  scale_y_continuous(limits = c(0, 65)) +
  scale_x_reverse(limits = rev(c(-100, 45000)))

# G/C skew
# RLSeq::plotEnrichment(rlr, )

# Genic Features
# table(sampleAnno$db)


# G4Q



### Figure 3/S3 ###
### RL Regions
dir.create("results/Figure_3", showWarnings = FALSE)

## Donut charts showing which samples went forward for analysis
dataNow <- rlsamples %>%
  filter(
    genome == "hg38",
    label == "POS",
    prediction == "POS",
    numPeaks > 5000,
    ip_type %in% c("S9.6", "dRNH")
  ) %>%
  group_by(ip_type) %>%
  tally()

ip_cols <- c("dRNH" = "#c386c4", "S9.6"="#82d0e8")

lt <- plot_ly(type = "pie") %>%
  add_pie(data=dataNow, labels = dataNow$ip_type, values = ~n, textinfo='label+value',
          marker = list(colors=ip_cols), 
          insidetextorientation='horizontal', hole=.6, rotation=250) %>%
  layout(showlegend = FALSE, title=list(text = "IP Type", x=0.15), margin = list(l = 100, r = 100, t=100, b=100),
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
save_image(plt, file = file.path(outpath, paste0("donut__", dat, ".svg")), format = "svg")


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
save_image(plt1, file = "results/Figure_3/dRNH__sample_donut.svg", format = "svg")
plt2 <- plot_ly(type = "pie") %>%
  add_pie(data=dataNow$S9.6, labels = dataNow$S9.6$Included, values = ~n, textinfo='label+value',
          marker = list(colors=ip_cols$S9.6), 
          insidetextorientation='horizontal', hole=.6) %>%
  layout(showlegend = FALSE, title=list(text = "Included", x=0.15), margin = list(l = 100, r = 100, t=100, b=100),
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
save_image(plt2, file = "results/Figure_3/S96__sample_donut.svg", format = "svg")


## R-loop venn diagram
rlregions <- rlregions_meta()
x <- rlregions %>%
  group_by(source) %>%
  {setNames(group_split(.), group_keys(.)[[1]])} %>% 
  lapply(function(x){
    pull(x, rlregion)
  })
x$dRNH <- c(x$dRNH, x$`dRNH S96`)
x$S96 <- c(x$S96, x$`dRNH S96`)
x$`dRNH S96` <- NULL
names(x) <- c("dRNH", "S9.6")
plt <- venn.diagram(x, fill = c("#e2a3e3", "#82d0e8"),
                    filename = NULL) %>%
  grid::grid.draw() %>%
  ggplotify::grid2grob() %>%
  ggplotify::as.ggplot()
ggsave(plt, filename = "results/Figure_3/venn_digram.svg")

## Get the significance via valr
drnh_cons <- read_tsv(file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_dRNH.narrowPeak"),
                      col_names = c("chrom", "start", "end", "name", "score", "strand", "signalVal", "pVal", "qVal", "peak"), 
                      show_col_types = FALSE, progress = FALSE)
s96_cons <- read_tsv(file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_S96.narrowPeak"),
                     col_names = c("chrom", "start", "end", "name", "score", "strand", "signalVal", "pVal", "qVal", "peak"), 
                     show_col_types = FALSE, progress = FALSE)
genome <- read_tsv("https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", 
                   col_names = c("chrom", "size"), show_col_types = FALSE, progress = FALSE)
print(valr::bed_fisher(drnh_cons, y = s96_cons, genome = genome))
valr::bed_intersect(x = drnh_cons, y = s96_cons) -> ddd
# Get sizes for venn diagram
length(which(! drnh_cons$name %in% ddd$name.x))
length(which(! s96_cons$name %in% ddd$name.y))
c(
  GenomicRanges::makeGRangesFromDataFrame(
    data.frame(chrom=ddd$chrom, start=ddd$start.x, end=ddd$end.x)
  ),
  GenomicRanges::makeGRangesFromDataFrame(
    data.frame(chrom=ddd$chrom, start=ddd$start.y, end=ddd$end.y)
  )
) %>% GenomicRanges::reduce()

# Density plot of peak sizes
pltdat <- bind_rows(list(
  tibble(ip_type = "dRNH", width = drnh_cons$end-drnh_cons$start),
  tibble(ip_type = "S9.6", width = s96_cons$end-s96_cons$start)
)) %>%
  filter(width < 80000) 
pltdat %>%
  group_by(ip_type) %>%
  summarise(
    density(width)$x[which.max(density(width)$y)]
  )
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
ggsave(plt, filename = "results/Figure_3/width_plot.svg")


### Figure 4/S4 ###
dir.create("results/Figure_4/", showWarnings = FALSE)

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

### ORganize RL Regions by row variance from the VST -- see what is the more variable

### PCA from RLRegion counts
cts <- RLHub::rlregions_counts()
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
signif(100*(pcdata$sdev / sum(pcdata$sdev)), 3)
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
pltpca
ggsave(plot = pltpca, filename = "results/Figure_4/pca_plot.svg", height = 5, width = 7)

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

# Volcano plot
volcplt <- restbl %>%
  EnhancedVolcano::EnhancedVolcano(
    lab = gsub(.$rlregion, pattern = "All_", replacement = ""), x = "log2FoldChange", y = "padj", 
    pCutoff = 1E-50, FCcutoff = 2,labSize = 4,
    title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none"
  ) +
  theme_bw(base_size = 14) +
  ggpubr::rremove("legend")
volcplt
ggsave(plot = volcplt, filename = "results/Figure_4/volcano.svg", height = 5, width = 6)
volcplt +
  theme_void() +
  ggpubr::rremove("legend") -> volcplt2
ggsave(plot = volcplt2, filename = "results/Figure_4/volcano.png", height = 5, width = 6)

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

enrichlist <- setNames(pbapply::pblapply(names(ol), function(groupnow) {
  genesNow <- ol[[groupnow]]
  response <- httr::POST(  # Post request to enrichr based on https://maayanlab.cloud/Enrichr/help#api&q=1
    url = 'https://maayanlab.cloud/Enrichr/addList', 
    body = list(
      'list' = paste0(genesNow, collapse = "\n"),
      'description' = groupnow
    )
  )
  response <- jsonlite::fromJSON(httr::content(response, as = "text"))  # Collect response
  permalink <- paste0("https://maayanlab.cloud/Enrichr/enrich?dataset=",  # Create permalink
                      response$shortId[1])
  return(permalink)
}), nm = names(ol))
enrichlist

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
ggsave(barplt, filename = "results/Figure_4/barplt.svg", height = 8, width = 8)


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
  ggsave(plot = plt, filename = paste0("results/Figure_4/pca_feature_plot_", pt,".svg"), height = 5, width = 7)
})


## Compare samples around genomic features
iptypeanno <- sampleAnno %>%
  left_join(rlsamples, by = c("experiment" = "rlsample")) %>%
  filter(! is.infinite(stat_fisher_rl), genome == "hg38",
         prediction == "POS", label == "POS",
         ip_type != "None", numPeaks > 5000,
         ! is.na(prediction), ! is.na(stat_fisher_rl)) %>%
  group_by(db, type) %>%
  mutate(stat_fisher_rl = (stat_fisher_rl+1) / (stat_fisher_shuf+1)) %>%
  select(c(ip_type, stat_fisher_rl)) %>%
  na.omit() %>%
  nest(data = c(ip_type, stat_fisher_rl)) %>%
  mutate(
    fit = map(data, ~ wilcox.test(stat_fisher_rl ~ ip_type, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(p.value, method = "bonferroni")
  ) %>%
  select(-data, -fit, -method, -alternative) %>%
  arrange(padj)
iptypeanno
pld <- RLSeq::plotEnrichment(rlr, pred_POS_only = TRUE, label_POS_only = TRUE, returnData = TRUE)
pld <- lapply(pld, function(x) {
  x %>%
    inner_join(rlsamples, by = c("experiment" = "rlsample")) %>%
    filter(ip_type %in% c("S9.6", "dRNH"), numPeaks > 5000)
})
hsplt <- pld$Encode_Histone
plot_multi_feature_iptype(features = unique(hsplt$type), pltdat = hsplt,
                          db = "Encode_Histone", lmts = c(-10, 15), 
                          factororder = unique(hsplt$type))
feats <- c("STAG1", "STAG2", "CTCF", "RAD21", "SMC3")
dbs <- c("Cohesin", "encodeTFBS")
hsplt <- bind_rows(as.list(pld[dbs])) 
plot_multi_feature_iptype(features = feats, pltdat = hsplt,
                          db = dbs, lmts = c(-10, 15), 
                          factororder = feats)
feats <- c("prom", "TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "polyA signal", "polyA site")
dbs <- c("Encode_CREs", "PolyA", "Transcript_Features")
hsplt <- bind_rows(as.list(pld[dbs])) 
plt <- plot_multi_feature_iptype(features = feats, pltdat = hsplt,
                          db = dbs, lmts = c(-10, 15), 
                          factororder = feats) + ggpubr::rotate_x_text(90)
plt +
  theme_classic(base_size = 14) -> plt
ggsave(plt, filename = "results/Figure_4/gene_feature_plot.svg", height = 6, width = 11)

## Do the same, but with regioneR
## TODO: Try this with summitted peaks
cons <- list(
  "dRNH" = drnh_cons,
  "S9.6" = s96_cons
)
if (! "annFull" %in% ls()) {
  annFull <- RLHub::annots_full_hg38()
}
feats <- c("enhD", "enhP", "prom", "TSS", "fiveUTR", "Exon", "Intron",
           "threeUTR", "TTS", "polyA_signal", "polyA_site",
           "STAG1", "STAG2", "CTCF", "RAD21", "SMC3", "MYC")
dbs <- c("Encode_CREs", "PolyA", "Transcript_Features",
         "Cohesin", "encodeTFBS")
keeps <- names(annFull)[grepl(names(annFull), pattern = paste0(dbs, collapse = "|")) &
                          grepl(names(annFull), pattern = paste0(feats, collapse = "|"))]
keepannots <- lapply(keeps, function(keep) {
  annFull[[keep]] %>% mutate(
    db = gsub(keep, pattern = "(.+)__(.+)", replacement = "\\1"),
    type = gsub(keep, pattern = "(.+)__(.+)", replacement = "\\2")
  )
}) %>% bind_rows() %>%
  filter(! (type == "CTCF" & db == "Encode_CREs")) %>%
  group_by(type) %>%
  {setNames(group_split(.), group_keys(.)[[1]])}
if (! file.exists("data/permTestRes.rds")) {
  permTestRes <- parallel::mclapply(keepannots, function(anno) {
    anno <- select(anno, -strand)
    annogr <- regioneR::toGRanges(as.data.frame(anno))
    annoRes <- setNames(lapply(cons, function(con) {
      con <- select(con, -strand)
      congr <- regioneR::toGRanges(as.data.frame(con))
      pt <- regioneR::permTest(A = congr, 
                               B = annogr, 
                               verbose = T, force.parallel = T, 
                               genome = "hg38", mask = RLSeq::genomeMasks$hg38.masked,
                               randomize.function = regioneR::circularRandomizeRegions,
                               evaluate.function = regioneR::numOverlaps)
      lz <- regioneR::localZScore(
        A = congr, B = annogr, pt, window = 5000, step = 50
      )
      return(list(
        "lz" = lz,
        "pt" = pt
      ))
    }), nm = names(cons))
    return(annoRes)
  }, mc.cores = length(keepannots))
  saveRDS(permTestRes, file = "data/permTestRes.rds")
} else {
  permTestRes <- readRDS("data/permTestRes.rds")
}

plts <- lapply(
  names(permTestRes),
  function(annonow) {
    rlfsReslst <- permTestRes[[annonow]]
    datNow <- lapply(
      names(rlfsReslst), function(nameNow) {
        rlfsResNow <- rlfsReslst[[nameNow]]
        # Obtain the plot data
        zs <- rlfsResNow$lz$`regioneR::numOverlaps`$shifted.z.scores
        pltdatNow <- dplyr::tibble(
          "zscore" = zs,
          "shift" = rlfsResNow$lz$`regioneR::numOverlaps`$shifts,
          "pval" = rlfsResNow$pt$`regioneR::numOverlaps`$pval,
          "ip_type" = nameNow
        )
      }
    ) %>% bind_rows() %>%
      group_by(ip_type) %>%
      # mutate(zscore = as.numeric(scale(zscore))) %>%
      # mutate(zscore = scales::rescale(zscore, to = c(0, 1))) %>%
      ungroup()
    s96p <- datNow$pval[datNow$ip_type == "S9.6"][1] %>% signif(5)
    drnhp <- datNow$pval[datNow$ip_type == "dRNH"][1] %>% signif(5)
    
    # Make plot
    pltbase <- ggplot2::ggplot(
      data = datNow,
      ggplot2::aes_string(
        y = "zscore",
        x = "shift",
        color = "ip_type"
      )
    ) +
      ggplot2::geom_line(size = 1) +
      ggplot2::ggtitle(annonow) +
      ggplot2::ylab("Z Score") +
      ggplot2::xlab("Shift") +
      ggprism::theme_prism(base_size = 14) 
  }
)
names(plts) <-  names(permTestRes)
plt <- ggpubr::ggarrange(plotlist = plts)


## Same but with deepTools
dir.create("data/annotation_beds", showWarnings = FALSE)
dir.create("data/annot_mats", showWarnings = FALSE)
pbapply::pbsapply(names(keepannots), function(keepannotsNow) {
  annotNow <- keepannots[[keepannotsNow]] %>% select(-strand)
  annobed <- paste0("data/annotation_beds/", annotNow$db[1], "__", annotNow$type[1], ".bed")
  annomat <- paste0("data/annot_mats/", annotNow$db[1], "__", annotNow$type[1], ".mat.gz")
  if (! file.exists(annobed)) {
    regioneR::toGRanges(as.data.frame(annotNow)) %>% 
      rtracklayer::export(con = annobed)
  }
  paste0(
    "computeMatrix reference-point -S ../RLBase-data/rlbase-data/rlregions/rlreg",
    "ions_S96.bw ../RLBase-data/rlbase-data/rlregions/rlregions_dRNH.bw ",
    "-R ", annobed, " -b 3000",
    " --referencePoint center -a 3000",
    " -o ", annomat, " -p 6 &"
  )
}) %>%
  write_lines(file = "data/annot_mats/anno_mat_cmd.sh")

## Run the script and then resume

# Plots
# plts <- parallel::mclapply(
#   names(keepannots), function(keepannotsNow) {
#     annots <- keepannots[[keepannotsNow]]
#     dtools_metaplot(db = annots$db[1], type = annots$type[1])
#   }, mc.cores = length(keepannots)
# )
# names(plts) <- names(keepannots)
# dir.create("results/Figure_4/plts")
# sapply(seq(plts), function(i) {
#   ggsave(plts[[i]], filename = paste0("results/Figure_4/plts/", names(plts)[i], ".png"))
# })
# plts <- ggpubr::ggarrange(plotlist = plts)


## Motifs
# sort -k5 -n -r ../RLBase-data/rlbase-data/rlregions/rlregions_S96.narrowPeak | head -n 15000 > data/S96_consensus_top.narrowPeak
# sort -k5 -n -r ../RLBase-data/rlbase-data/rlregions/rlregions_dRNH.narrowPeak | head -n 15000 > data/dRNH_consensus_top.narrowPeak
# fastaFromBed -fi ~/.rlpipes_genomes/hg38/hg38.fa -bed data/dRNH_consensus_top.narrowPeak > data/dRNH_consensus.fasta
# fastaFromBed -fi ~/.rlpipes_genomes/hg38/hg38.fa -bed data/S96_consensus_top.narrowPeak > data/S96_consensus.fasta
# meme-chip -oc data/dRNH_meme_res_dif -db data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -meme-nmotifs 0 -neg data/S96_consensus.fasta data/dRNH_consensus.fasta
# meme-chip -oc data/S96_meme_res_dif -db data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -meme-nmotifs 0 -neg data/dRNH_consensus.fasta data/S96_consensus.fasta

## Genes
genes <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, 
  keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86), 
  columns = c("SEQNAME", "GENESEQSTART", "GENESEQEND", "SYMBOL")
) %>% dplyr::rename(chrom = SEQNAME, start = GENESEQSTART, end = GENESEQEND) %>%
  filter(! grepl(GENEID, pattern = "LRG.+")) %>%
  select(-GENEID) %>% mutate(chrom = paste0("chr", chrom)) %>%
  distinct(SYMBOL, .keep_all = TRUE)
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

# Annotate TTS, TSS, and Gene Body
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
  scale_fill_manual(values = setNames(
    rev(c("#b36262", "#c4cc6a", "#d17bdb", "#4bc6d6", "#83d647", "#9494e3", "#7d7d7d")),
    nm = rev(c("TSS", "fiveUTR", "Exon", "Intron", "threeUTR", "TTS", "Intergenic"))
  ), guide = guide_legend(title = "Feature", reverse = TRUE))
plt

# Plot gene distribution
# GO
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

ggplot(pltdats, aes(x = rank, y = pct_cons, color = source, label=label)) +
  geom_hline(yintercept = 50, color = "grey", linetype = "dashed", alpha = .85) +
  geom_point() +
  ggrepel::geom_label_repel(max.overlaps = Inf, color = "black", 
                            force = 50, nudge_x = -500, nudge_y = -.05) +
  facet_wrap(~ source, scales = "free_x") +
  theme_bw(base_size = 14) +
  ylab("RL region conservation (% samples)") +
  xlab("RL region rank")

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
  slice_max(n = 1000, order_by = pct_cons) %>%
  mutate(SYMBOL = strsplit(SYMBOL, split = "\n")) %>%
  unnest(SYMBOL) %>%
  filter(SYMBOL %in% ribogenes)
complst <- data.frame(
  "S9.6" = ribogenes %in% unique(toplt$SYMBOL[toplt$source == "S9.6"]),
  "dRNH" = ribogenes %in% unique(toplt$SYMBOL[toplt$source == "dRNH"]),
  "Ribo genes" = rep(TRUE, length(ribogenes))
) %>% as.matrix()

fit <- eulerr::euler(complst)
plot(fit,
     quantities = TRUE,
     lty = 1,
     labels = list(font = 4))

rnks <- pltdats$rank[pltdats$source == "S9.6"]
x <- pltdats$rank[pltdats$source == "S9.6" & pltdats$is_ribo]
rs <- ks.test(x = x, y = sample(rnks, length(x)))
  
rnks <- pltdats$rank[pltdats$source == "dRNH"]
x <- pltdats$rank[pltdats$source == "dRNH" & pltdats$is_ribo]
rs2 <- ks.test(x = x, y = sample(rnks, length(x)))

# Show this on genome browser
pltdats %>% 
  mutate(SYMBOL = strsplit(SYMBOL, split = "\n")) %>%
  unnest(SYMBOL) %>%
  filter(SYMBOL %in% ribogenes) %>%
  select(SYMBOL, pct_cons, source) %>%
  pivot_wider(id_cols = SYMBOL, values_from = pct_cons, names_from = source, values_fill = list(0)) %>%
  unnest(dRNH) %>%
  unnest(S9.6) %>% mutate(diff = S9.6 - dRNH) -> dd

# Run permTest to show this
anno <- genes %>% filter(SYMBOL %in% ribogenes)
annogr <- regioneR::toGRanges(as.data.frame(anno))
rtracklayer::export(annogr, con = "data/ribogenes.bed")
metap <- read_tsv("data/ribogenes_matrix.mat.gz", skip = 1, col_names = FALSE)
metadat <- read_lines("data/ribogenes_matrix.mat.gz", n_max = 1)
metadat <- gsub(metadat, pattern = "@", replacement = "") %>% jsonlite::fromJSON()
mat <- metap %>%
  select(-c(1:6)) 
mat %>%
  colMeans(na.rm = TRUE) -> cm
divs <- setNames(c(42, 184), nm = c("rlregions_dRNH", "rlregions_S96"))
pltdats <- tibble("score" = cm,
                  "loc" = rep(seq(length(cm) / 2), 2),
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
tss <- (3/11)*length(cm)/2
tts <- tss + (5/11)*length(cm)/2
plt <- pltdats %>%
  ggplot(aes(x = loc, y = score, color = label, fill=label)) +
  geom_vline(xintercept = tss,  color="grey") +
  geom_vline(xintercept = tts,  color="grey") +
  geom_line() +
  theme_classic(base_size = 14) +
  ylab("Density") +
  xlab("Gene position (bp)") +
  scale_x_continuous(
    breaks = c(1, tss, tts, length(cm)/2),
    labels = c("-3000", "TSS", "TTS", "+3000"),
    expand = c(0,0)
  ) +
  ggpubr::rremove("legend") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(mcm, mxcm),
                     breaks = mcm + c(0*mcdiff, .25*mcdiff, .5*mcdiff, .75*mcdiff, 1*mcdiff),
                     labels = c(0, .25, .5, .75, 1)) 
plt
ggsave(plt, filename = "results/Figure_4/gene_position_rlregions.svg", height = 4.35, width = 6.0)

annoRes <- setNames(pbapply::pblapply(cons, function(con) {
  con <- select(con, -strand)
  congr <- regioneR::toGRanges(as.data.frame(con))
  pt <- regioneR::permTest(A = congr, 
                           B = annogr, 
                           verbose = T, force.parallel = T, 
                           genome = "hg38", mask = RLSeq::genomeMasks$hg38.masked,
                           randomize.function = regioneR::circularRandomizeRegions,
                           evaluate.function = regioneR::numOverlaps)
  lz <- regioneR::localZScore(
    A = congr, B = annogr, pt, window = 5000, step = 50
  )
  return(list(
    "lz" = lz,
    "pt" = pt
  ))
}), nm = names(cons))


# For each gene, do a comparison plot
pltdats2 <- pltdats %>%
  top_frac(n = .10, wt = pct_cons) %>%
  mutate(SYMBOL = strsplit(SYMBOL, split = "\n")) %>%
  unnest(SYMBOL) 
complst <- pltdats2 %>%
  {setNames(group_split(.), nm = group_keys(.)[[1]])} %>%
  lapply(function(x) unique(pull(x, "SYMBOL")))
plt <- venn.diagram(complst, fill = c("#e2a3e3", "#82d0e8"),
                    filename = NULL) %>%
  grid::grid.draw() %>%
  ggplotify::grid2grob() %>%
  ggplotify::as.ggplot()

ol <- calculate.overlap(complst)
ol <- list(
  "dRNH-only" = complst$dRNH[! complst$dRNH %in% complst$S9.6],
  "S9.6-only" = complst$S9.6[! complst$S9.6 %in% complst$dRNH],
  "Shared" = complst$dRNH[complst$dRNH %in% complst$S9.6]
)
enrichlist <- setNames(pbapply::pblapply(names(ol), function(groupnow) {
  genesNow <- ol[[groupnow]]
  response <- httr::POST(  # Post request to enrichr based on https://maayanlab.cloud/Enrichr/help#api&q=1
    url = 'https://maayanlab.cloud/Enrichr/addList', 
    body = list(
      'list' = paste0(genesNow, collapse = "\n"),
      'description' = groupnow
    )
  )
  response <- jsonlite::fromJSON(httr::content(response, as = "text"))  # Collect response
  permalink <- paste0("https://maayanlab.cloud/Enrichr/enrich?dataset=",  # Create permalink
                      response$shortId[1])
  return(permalink)
}), nm = names(ol))

# BReak by location
pltlst <- pltdats %>%
  mutate(SYMBOL = strsplit(SYMBOL, split = "\n")) %>%
  unnest(SYMBOL)  %>%
  mutate(type2 = ifelse(type %in% c("TTS", "TSS", "Intergenic"), type, "Gene Body")) %>%
  filter(type != "Intergenic") %>%
  mutate(source2 = paste0(source, type2)) %>%
  group_by(source2) %>%
  slice_max(n = 1000, order_by = pct_cons) %>%
  group_by(type2) %>%
  {setNames(group_split(.), nm = group_keys(.)[[1]])}
plts <- pltlst %>%
  lapply(function(x) {
    complst <- x %>%
      group_by(source) %>%
      {setNames(group_split(.), nm = group_keys(.)[[1]])} %>%
      lapply(function(x) unique(pull(x, "SYMBOL")))
    venn.diagram(complst, fill = c("#e2a3e3", "#82d0e8"),
                 filename = NULL) %>%
      grid::grid.draw() %>%
      ggplotify::grid2grob() %>%
      ggplotify::as.ggplot()
  })
ggpubr::ggarrange(plotlist = plts)
# enrich
pltlst %>% lapply(function(x) {
  x %>% group_by(source) %>% {setNames(group_split(.), nm = group_keys(.)[[1]])}
}) -> toEnrich
enrichlist <- setNames(pbapply::pblapply(names(toEnrich), function(groupnow) {
  sources <- toEnrich[[groupnow]]
  sapply(sources, function(x) {
    groupfinal <- paste0(x$source[1], "__", groupnow)
    genesNow <- unique(x$SYMBOL)
    response <- httr::POST(  # Post request to enrichr based on https://maayanlab.cloud/Enrichr/help#api&q=1
      url = 'https://maayanlab.cloud/Enrichr/addList', 
      body = list(
        'list' = paste0(genesNow, collapse = "\n"),
        'description' = groupfinal
      )
    )
    response <- jsonlite::fromJSON(httr::content(response, as = "text"))  # Collect response
    paste0("https://maayanlab.cloud/Enrichr/enrich?dataset=",  # Create permalink
           response$shortId[1])
  })
  
}), nm = names(toEnrich))




# Redo meme chip with summitting
pltdats2 %>%
  select(chrom, start, end, name, source) %>%
  {setNames(group_split(.), nm = group_keys(.)[[1]])} %>%
  lapply(function(x) {
    src <- x$source[1]
    res <- as.data.frame(x) %>%
      GenomicRanges::makeGRangesFromDataFrame() %>%
      GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
      GenomicRanges::trim() %>%
      unique()
    res %>%
      rtracklayer::export(con = paste0("data/", src, "_consensus_top.bed"))
    return(res)
  })




### With Case vs Control vs RNH/nonRNH confusion matrix
rmapfftsmall %>%
  mutate(Condition = ifelse(is_rnh_like, "RNH-like", "Normal"),
         Prediction = ifelse(prediction == "case", "Normal", "RNH-like")) %>%
  group_by(Condition, Prediction) %>%
  tally() %>%
  pivot_wider(id_cols = c(Prediction), names_from = Condition, values_from = n, values_fill = 0) %>%
  # column_to_rownames(var = "Predicted") %>%
  kable(format = "latex", booktabs=T, 
        caption = "Predicted vs Labeled Condition", escape = FALSE, align = "c") %>%
  kableExtra::add_header_above(c(" " = 1, "Label" = 2))

toTbl <- rmapfftsmall %>%
  mutate(Condition = ifelse(is_rnh_like, "RNH-like", "Normal"),
         Prediction = ifelse(prediction == "case", "Normal", "RNH-like")) %>%
  group_by(Condition, Prediction, mode) %>%
  tally() %>%
  pivot_wider(id_cols = c(mode),
              names_from = c(Prediction, Condition),
              values_from = n, names_sep = NAME_SEP,
              names_sort = TRUE, values_fill = 0) %>%
  column_to_rownames("mode")
cnms1 <- gsub(colnames(toTbl), 
              pattern = paste0('(.+)', NAME_SEP, "(.+)"), replacement = "\\1")
cnms2 <- gsub(colnames(toTbl), 
              pattern = paste0('(.+)', NAME_SEP, "(.+)"), replacement = "\\2")
toMerge <- match(unique(cnms1), cnms1)
mergeVec <- sapply(seq(toMerge), function(i) {
  if (i == length(seq(toMerge))) {
    ret <- toMerge[i]:length(cnms1)
    length(ret) 
  } else {
    ret <- toMerge[i]:toMerge[i+1]
    length(ret[-length(ret)])
  }
})
names(mergeVec) <- unique(cnms1)
mergeVec <- c("Predicted" = 1, mergeVec)

cbind(colSums(toTbl),
      rep(" ", length(colSums(toTbl)))) %>%
  t() %>%
  cbind(data.frame(names = c("Sum", "Mode"))) %>%
  mutate(names = cell_spec(
    names, format = "latex", bold = names == "Mode"
  )) %>%
  column_to_rownames(var = "names") %>%
  rbind(toTbl) %>%
  rownames_to_column("Label") %>%
  kable(toTbl, 
        format = "latex", 
        caption = "Number of samples by Mode, Condition, Prediction",
        escape = FALSE,
        align = "r",
        booktabs = T,
        col.names = c("Labeled",cnms2 )) %>%
  add_header_above(mergeVec)

### Compare RLFS pval between modes
rlfs_sig <- dataLst %>%
  pluck("sample_quality_characteristics") %>%
  filter(char_type == "rlfs_pval") %>%
  mutate(RLFS_sig = ifelse(value > 1.3, "p < .05", "n.s.")) 
rmapfft_rlfs <- left_join(rmapfft_peaks, rlfs_sig)

toTbl <- rmapfft_rlfs %>%
  group_by(mode, RLFS_sig) %>%
  tally() %>%
  pivot_wider(id_cols = c(mode), names_from = RLFS_sig, values_from = n, values_fill = 0) %>%
  column_to_rownames("mode")

cbind(colSums(toTbl),
      rep(" ", length(colSums(toTbl)))) %>%
  t() %>%
  cbind(data.frame(names = c("Sum", "Mode"))) %>%
  mutate(names = cell_spec(
    names, format = "latex", bold = names == "Mode"
  )) %>%
  column_to_rownames(var = "names") %>%
  rbind(toTbl) %>%
  kable(format = "latex", booktabs=T, 
        caption = "RLFS-enrichment testing results by mode", escape = FALSE, align = "c") %>%
  kableExtra::add_header_above(c(" " = 1, "RLFS-enrichment" = 2))


### Compare RLFS pval between modes


### Compare # peaks between modes


### Feature enrichment analysis

## Build the table of annotation types with the 
## Significance of enrichment (wilcox) within R-loop mapping samples
rlfeatres <- sampleAnno %>%
  left_join(rlsamples, by = c("experiment" = "rlsample")) %>%
  filter(! is.infinite(stat_fisher_rl), genome == "hg38",
         ! is.na(prediction), ! is.na(stat_fisher_rl)) %>%
  group_by(db, type) %>%
  mutate(stat_fisher_rl = (stat_fisher_rl+1) / (stat_fisher_shuf+1)) %>%
  select(c(prediction, stat_fisher_rl)) %>%
  na.omit() %>%
  nest(data = c(prediction, stat_fisher_rl)) %>%
  mutate(
    fit = map(data, ~ wilcox.test(stat_fisher_rl ~ prediction, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(p.value, method = "bonferroni")
  ) %>%
  select(-data, -fit, -method, -alternative) %>%
  arrange(padj)
write_csv(rlfeatres, file.path("results", "tables", "RLBase_anno_enrichment_summarized.csv"))

## Show the strength of these results with the top hit
rlr <- RLSeq::RLRangesFromRLBase("SRX1070676")
rlr <- RLSeq::featureEnrich(rlr, annotype = "full")
pltdat <- RLSeq::plotEnrichment(rlr, returnData = TRUE, 
                                pred_POS_only = FALSE,
                                rlbaseRes = sampleAnno, 
                                rlsamples = rlsamples,
                                label_POS_only = FALSE) 
pltdat$G4Qpred <- pltdat$G4Qpred %>%
  filter(stat_fisher_rl != -10, stat_fisher_rl != 15)
plot_feature(feature = "tl:4-5 nl:1-2 gn:1", db = "G4Qpred", 
             pltdat = pltdat,
             lmts = c(-6, 10),
             yb = 8.5) +
  ggtitle("Feature Enrichment Results in RLBase samples",
          subtitle = "G4 Quadruplex Predicted sites (4-5:1-2:1)")


## Compare with RLBPs
rlbpGeneList <- rlbps$geneName[rlbps$geneName %in% humanGenes$SYMBOL & rlbps$combinedScore > 0]
rlfeatres %>%
  filter(
    type %in% humanGenes$SYMBOL
  ) %>%
  mutate(
    RLBP = type %in% {{ rlbpGeneList }},
    padj = -log10(padj),
    statisticrank = rank(-statistic),
    label = case_when(
      RLBP & padj > 30 ~ type, TRUE ~ ""
    )
  ) %>%
  arrange(RLBP) %>%
  ggplot(aes(x = statistic, y = padj, color = RLBP, label=label)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 10, force = 5,
                           nudge_y = -1, max.time = 5, max.iter = 20000, nudge_x = 100) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "#dedede"))  +
  geom_rug(sides = "b") +
  xlab("Wilcox Statistic") + ylab("P Adjusted Value (-log10)") +
  ggtitle("R-loop binding factors vs known RLBPs") +
  theme_bw(base_size = 14) +
  scale_y_continuous(limits = c(0, 65)) +
  scale_x_reverse(limits = rev(c(-100, 45000)))

## Compare with RNAi
RNAi <- read_csv("data/rnai_2013_results.csv")
rnaiGeneList <- RNAi$geneName[RNAi$geneName %in% humanGenes$SYMBOL & RNAi$zScore < -2]
rlfeatres %>%
  filter(
    type %in% humanGenes$SYMBOL
  ) %>%
  mutate(
    RNAi = type %in% {{ rnaiGeneList }},
    padj = -log10(padj),
    label = case_when(
      RNAi & padj > 25 ~ type, TRUE ~ ""
    )
  ) %>%
  arrange(RNAi) %>%
  ggplot(aes(x = statistic, y = padj, color = RNAi, label=label)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 30, force = 2) +
  scale_color_manual(values = c("TRUE" = "#284880", "FALSE" = "#dedede"))  +
  scale_x_reverse() + geom_rug(sides = "b") +
  xlab("Wilcox Statistic") + ylab("P Adjusted Value (-log10)") +
  ggtitle("R-loop binding factors vs RNAi hits") +
  theme_bw(base_size = 14)

## Show the difference between the DNA binding and RNA binding
dups <- rlfeatres$type[duplicated(rlfeatres$type)]
rlfeatres %>%
  filter(db %in% c("RBP_ChIP", "encodeTFBS")) -> chips
rlfeatres %>%
  filter(db %in% "RBP_eCLiP", type %in% chips$type) -> keep
rlfeatres %>%
  filter(db %in% c("RBP_ChIP", "RBP_eCLiP", "encodeTFBS"),
         type %in% keep$type) %>%
  mutate(moeity = case_when(db != "RBP_eCLiP" ~ "DNA (ChIP)", TRUE ~ "RNA (eCLiP)"),
         score = -log10(padj)) %>%
  group_by(type, moeity) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  pivot_wider(id_cols = "type", names_from = "moeity", values_from = "score") %>%
  mutate(differ = `RNA (eCLiP)` - `DNA (ChIP)`) %>%
  pivot_longer(cols = 2:3) %>%
  mutate(padj = abs(value)) %>% 
  arrange(differ) %>%
  ungroup() %>% 
  mutate(type = factor(type, levels = unique(type))) %>%
  ggplot(aes(x = type, y = padj, fill = name)) +
  geom_col(width = .7, position = position_dodge(.7)) +
  theme_prism(base_size = 14) + 
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0)) +
  ggpubr::rotate_x_text(30, vjust = 1) +
  ggtitle("R-loop associated factors with eCLiP and CHIP enrichment") +
  xlab("Factor") + ylab("Enrichment Padj value (-log10)") +
  scale_fill_manual(values = c("goldenrod", "navy"))

## Show the enrichment results for the RL Regions split by dRNH and S96
rlrAnno %>%
  filter(! is.infinite(stat_fisher_rl), 
         ! is.na(opt), ! is.na(stat_fisher_rl)) %>%
  filter(opt != "All") %>%
  relocate(opt, .before = db) %>%
  group_by(db, type) %>%
  filter(all(c("dRNH", "S96") %in% opt)) %>%
  mutate(stat_fisher_rl2 = (stat_fisher_rl+1) / (stat_fisher_shuf+1)) %>%
  select(c(opt, stat_fisher_rl2)) %>%
  na.omit() %>%
  ungroup() %>%
  pivot_wider(id_cols = c("db", "type"), names_from = opt, values_from = "stat_fisher_rl2") %>%
  mutate(dRNH_scale = scale(dRNH)[,1],
         S96_scale = scale(S96)[,1]) %>%
  group_by(db, type) %>%
  mutate(differ = sd(c(dRNH, S96)),
         differ_scale = sd(c(dRNH_scale, S96_scale))) %>%
  arrange(desc(differ_scale)) %>% 
  ungroup() %>%
  mutate(top = row_number() <= 100,
         label= ifelse(top, type, "")) %>%
  filter(db %in% c("Encode_Histone", "encodeTFBS", "G4Qpred__G4Pred", "RBP_ChIP", "RBP_eCLiP")) %>%
  ggplot(aes(x = S96_scale, y = dRNH_scale, color=db, label=label)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  ggrepel::geom_text_repel() +
  scale_y_continuous(limits = c(-1, 5)) +
  scale_x_continuous(limits = c(-1, 5)) +
  # scale_y_log10(limits = c(0.01, 10)) +
  # scale_x_log10(limits = c(0.01, 10)) +
  ylab("dRNH scaled enrichment") +
  xlab("S9.6 scaled enrichment") + 
  labs(title = "Feature Enrichment of dRNH vs S9.6 consensus RLoops")

