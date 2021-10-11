### Main script for generating figures (eventually should be notebook) ###

library(tidyverse)
library(kableExtra)
library(tidymodels)
library(VennDiagram)
library(ggprism)

# Load data from RLHub
sampleAnno <- RLHub::feat_enrich_samples()
rlrAnno <- RLHub::feat_enrich_rlregions()
rlsamples <- RLHub::rlbase_samples()
rlbps <- RLHub::rlbps()

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

### RLBase data

## Get the summary of the datasets ##

# Breakdown of dataset categorical variables
toTbl <- rlsamples %>%
  group_by(mode, genome) %>%
  tally() %>%
  pivot_wider(id_cols = genome, names_from = mode, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "genome") %>%
  t() 
cbind(colSums(toTbl),
      rep(" ", length(colSums(toTbl)))) %>%
  t() %>%
  cbind(data.frame(names = c("Sum", "Mode"))) %>%
  mutate(names = cell_spec(
    names, format = "html", bold = names == "Mode"
  )) %>%
  column_to_rownames("names") %>%
  rbind(toTbl) %>%
  kable(format = "html", booktabs=T, 
        caption = "Number of samples by Mode, Genome",
        escape = FALSE, align = "r") %>%
  kableExtra::kable_styling()


NAME_SEP <- "__xXx__"
toTbl <- rlsamples %>%
  group_by(mode, genome, label) %>%
  tally() %>%
  pivot_wider(id_cols = c(mode),
              names_from = c(genome, label),
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
mergeVec <- c(" " = 1, mergeVec)

cbind(colSums(toTbl),
      rep(" ", length(colSums(toTbl)))) %>%
  t() %>%
  cbind(data.frame(names = c("Sum", "Mode"))) %>%
  mutate(names = cell_spec(
    names, format = "html", bold = names == "Mode"
  )) %>%
  column_to_rownames("names") %>%
  rbind(toTbl) %>%
  kable(toTbl, 
        format = "html", 
        caption = "Number of samples by Mode, Condition, Genome",
        escape = FALSE,
        align = "r",
        booktabs = T,
        col.names = cnms2) %>%
  add_header_above(mergeVec) %>%
  kableExtra::kable_styling()

## Scatter pie charts
# Mode
forPie <- rlsamples %>%
  mutate(mode = case_when(
    mode %in% RLSeq:::auxdata$mode_cols$mode ~ mode,
    TRUE ~ "misc"
  )) %>%
  group_by(mode) %>%
  tally() %>%
  left_join(RLSeq:::auxdata$mode_cols, by = "mode") %>%
  mutate(label = paste0(
    mode, "\n", round(100*n/sum(n), 1), "%\n(", n, ")"
  )) %>%
  sample_frac()
pie(forPie$n, labels = forPie$label, 
    clockwise = F, init.angle = 150, cex = .77,  
    col = forPie$col)

# Genome
forPie <- rlsamples %>%
  group_by(genome) %>%
  tally() %>%
  mutate(label = paste0(
    genome, "\n", round(100*n/sum(n), 1), "%\n(", n, ")"
  )) %>%
  sample_frac()
pie(forPie$n, labels = forPie$label, 
    clockwise = F, init.angle = 150, cex = .77)


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
ComplexHeatmap::pheatmap(heatdata, color = colors1, main = "Percentage samples labels by Mode",
                         fontsize = 13, cluster_rows = FALSE, cluster_cols = FALSE,
                         display_numbers = TRUE, fontsize_number = 10.5,
                         annotation_col = annocol,
                         show_colnames = FALSE, name = "Representation (%)",
                         annotation_colors = hm$cat_cols[c("label", "prediction")],
                         annotation_row = annorow[rownames(heatdata),,drop=FALSE],
                         number_color = "black", angle_col = "0")


# Num peaks across conditions 
rlsamples %>%
  filter(prediction == "POS") %>%
  mutate(numPeaks = numPeaks / 1000) %>%
  ggplot(mapping = aes(x = mode, fill = mode,
                       y = numPeaks)) +
  geom_boxplot(color = "black", width = .7) +
  geom_jitter() +
  ylab("Peaks Called (thousands)") +
  xlab(NULL) +
  labs(title = "Peaks called per mode", subtitle = "Without 'NEG' predicted samples") +
  scale_fill_manual(values = setNames(RLSeq:::auxdata$mode_cols$col, nm = RLSeq:::auxdata$mode_cols$mode)) +
  guides(y = "prism_offset_minor") + 
  theme_prism(base_size = 16) + 
  ggpubr::rremove("legend") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1)) +
  ggpubr::stat_compare_means(size = 6)


# Num of peaks across condition + RNH status
DODGE_POS <- .8
df_p_val <- rstatix::group_by(rlsamples, mode) %>%
  filter(! is.na(prediction), length(which(prediction == "POS")) > 2,
         length(which(prediction == "NEG")) > 2) %>%
  rstatix::wilcox_test(numPeaks ~ prediction, p.adjust.method = "fdr") %>%
  rstatix::add_significance(p.col = "p") %>% 
  rstatix::add_xy_position(dodge = DODGE_POS, x = "mode")

rlsamples %>%
  right_join(df_p_val, by = "mode") %>%
  ggplot(mapping = aes(x = mode, 
                       y = numPeaks, fill = label)) +
  geom_boxplot(color = "black", width = .7,
               position = position_dodge(DODGE_POS)) +
  geom_jitter(position = position_jitterdodge(.2)) +
  ylab("Peaks Called (thousands)") +
  xlab(NULL) +
  labs(title = "Peaks called per mode", subtitle = "+/- RNH-like samples") +
  scale_fill_prism("prism_light", name = "RNaseH1") +
  guides(y = "prism_offset_minor") + 
  theme_prism(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1)) +
  add_pvalue(df_p_val,  
             xmin = "xmin", 
             xmax = "xmax",
             inherit.aes = F)


### PCA from RLRegion counts
cts <- RLHub::rlregions_counts()
ctsPos <- cts[,! is.na(cts$prediction) & cts$prediction == "POS"]
ctsmat <- ctsPos@assays@data$cts
vstmat <- DESeq2::vst(ctsmat)
pcdata <- prcomp(t(vstmat))
pcsmall <- pcdata$x[,1:7]
tsdat <- Rtsne(pcsmall, dims = 2, pca = FALSE, theta = 0)
tdat <- as.data.frame(tsdat$Y)
rownames(tdat) <- rownames(pcdata$x)  
colnames(tdat) <- c("TSNE_1", "TSNE_2")
pltdat <- tdat %>% 
  rownames_to_column(var = "name") %>%
  inner_join(rlsamples, by = c("name" = "rlsample")) %>%
  mutate(mode = case_when(
    mode %in% RLSeq:::auxdata$mode_cols$mode ~ mode,
    TRUE ~ "misc"
  ))
pltdat %>%
  {
    ggplot(., aes(x = TSNE_1, y = TSNE_2, color = mode, group = tissue, label = name)) +
      geom_point() +
      scale_color_manual(values = setNames(RLSeq:::auxdata$mode_cols$col,
                                           nm = RLSeq:::auxdata$mode_cols$mode)) +
      theme_bw(base_size = 14) +
      ggtitle("TSNE plot of RLBase samples")
  } 

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
plot_feature <- function(feature, db, pltdat, lmts, yb) {
  
  # Define limits
  splitby <- "prediction"
  cols <- RLSeq:::auxdata[[paste0(tolower(splitby[1]), "_cols")]]
  colvec <- cols %>% dplyr::pull(.data$col)
  names(colvec) <- as.data.frame(cols)[, splitby[1]]
  yval <- "stat_fisher_rl"
  pltdat %>%
    as.list() %>%
    bind_rows() %>%
    filter(type == {{ feature }}, 
           db == {{ db }}) %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = "prediction",
        y = yval,
        label = "experiment",
        fill = if (splitby[1] != "none") splitby[1],
        color = if (splitby[1] != "none") splitby[1]
      )
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_violin(
      color = "black",
      width = 1,
      trim = FALSE,
      alpha = .5,
      na.rm = TRUE
    ) +
    ggplot2::geom_boxplot(
      width = .5,
      color = "black",
      alpha = 1
    ) +
    ggplot2::scale_fill_manual(
      values = colvec, drop = TRUE
    ) + 
    ggplot2::ylab(
      "Fisher's test odds ratio (log2)"
    ) +
    ggplot2::xlab("RLSeq Prediction") +
    ggplot2::labs(
      title = paste0(feature, " - ", db),
      subtitle = "Feature enrichment results"
    ) +
    ggpubr::stat_compare_means(comparisons = list(c("POS", "NEG")),
                               label.y = yb,
                               label = "p.signif", size=6) +
    theme_prism(base_size = 14) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 14),
      plot.caption = ggplot2::element_text(size = 11)
    ) +
    ggpubr::rremove("legend") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(size = 0, stroke = 0)
      )
    ) +
    scale_y_continuous(limits = lmts)
}
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

