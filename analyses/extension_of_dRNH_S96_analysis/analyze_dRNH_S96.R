## This is an extension of the analysis from the original paper 
## Purpose is to characterize S96 and dRNH peaks and dig deeper into them

library(tidyverse)

#### MAGIC ####

ndrnh <- 56  # Number of dRNH passing QC filters
ns96 <- 191  # Number of S9.6 passing QC filters

###############


drnh_cons <- read_tsv(file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_dRNH.narrowPeak"),
                      col_names = c("chrom", "start", "end", "name", "score", "strand", "signalVal", "pVal", "qVal", "peak"), 
                      show_col_types = FALSE, progress = FALSE)
s96_cons <- read_tsv(file.path(RLSeq:::RLBASE_URL, "misc", "rlregions_S96.narrowPeak"),
                     col_names = c("chrom", "start", "end", "name", "score", "strand", "signalVal", "pVal", "qVal", "peak"), 
                     show_col_types = FALSE, progress = FALSE)
genome <- read_tsv("https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", 
                   col_names = c("chrom", "size"), show_col_types = FALSE, progress = FALSE)


## Get the dRNH-only and S96-only along with shared
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
  filter(width < 80000) 
plt <- ggplot(pltdat, aes(x = width, color = group, fill = group)) +
  geom_density(alpha=.3, adjust=2) +
  scale_x_log10(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic(base_size = 18) +
  ylab("Frequency") + xlab("Peak width (bp)") +
  guides(fill=guide_legend(title=NULL), color=guide_legend(title=NULL))



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

