#' Makes the "RLoopSignal" table for RMapDB
#'
#' @importFrom magrittr %>%
#' @import rlang
makeRLoopSignal <- function() {
  # TODO: Should we consider a separate category for promoter overlap?
  ChIPpeakAnno::findOverlapsOfPeaks(peakLst, ignore.strand = FALSE) %>%
    purrr::pluck("overlappingPeaks") %>%
    purrr::pluck(1) %>%
    dplyr::select(
      !! names(peakLst)[1] := peaks1,
      !! names(peakLst)[2] := peaks2,
    )
}

#' Helper function for converting CSV to GR
csvToGR <- function(csv) {
  readr::read_csv(csv) %>%
    dplyr::mutate(
      seqnames = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\1"),
      start = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\2"),
      end = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\3"),
      strand = gsub(location, pattern = "(.+):(.+)\\-(.+):(.+)", replacement = "\\4"),
      # from https://stackoverflow.com/questions/45146688/execute-dplyr-operation-only-if-column-exists
      score = if ("confidence_level" %in% colnames(.)) confidence_level else 0
    ) %>%
    dplyr::filter(! grepl(pattern = "LRG.+", x = id)) %>%
    dplyr::select(
      seqnames,
      start,
      end,
      name=id,
      score,
      strand
    ) %>%
    dplyr::distinct(seqnames, start, end, .keep_all = TRUE) %>%
    as.data.frame() %>%
    # Resolves issue of having NAs in the start/end column
    dplyr::mutate(
      start = as.numeric(as.character(start)),
      end = as.numeric(as.character(end))
    ) %>%
    na.omit() %>%
    ChIPpeakAnno::toGRanges() 
}

if (sys.nframe() == 0L) {
  source("analyses/Prepare-RMapDB-Tables/GeneRLoopOverlap_GenFeatRLoopOverlap/gene_rl_overlap__gf_rl_overlap.R")
  
  # TODO: Needs to become part of the snake pipe
  RLOOPS <- "analyses/Prepare-RMapDB-Tables/rloops.csv.xz"
  RMAPSAMPS <- "analyses/Prepare-RMapDB-Tables/rmap_samples.csv.xz"
  RLOOP_SIGNAL_MAT <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloop_signal.tsv.xz"
  
  # Need to create BED file for deeptools to use
  RLOOPS_BED <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloops.bed"
  RLOOPS %>% csvToGR() %>% rtracklayer::export.bed(con = RLOOPS_BED)
  system(paste0("xz -f ", RLOOPS_BED))
  
  # Get RLoops
  rloops <- RLOOPS %>% csvToGR() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "names") %>%
    dplyr::mutate(loc = paste0(seqnames, "_", start, "_", end)) %>%
    dplyr::select(-c("seqnames", "start", "end", "width", "strand", "score"))
  
  # Read in the TSV
  read_tsv(RLOOP_SIGNAL_MAT) %>% 
    dplyr::mutate(loc = paste0(`#'chr'`, "_", `'start'` + 1, "_", `'end'`)) %>%
    dplyr::select(-c(`#'chr'`, `'start'`,`'end'`)) %>%
    right_join(rloops, by = "loc") %>%
    dplyr::select(-loc) %>%
    dplyr::rename(rloop_id=names) %>%
    tidyr::pivot_longer(cols = ! rloop_id,
                        names_to = "rmap_sample_id",
                        values_to = "norm_counts") %>%
    dplyr::mutate(rmap_sample_id = gsub(rmap_sample_id,
                                        pattern = "'(SRX[0-9]+)_.+",
                                        replacement = "\\1"))
  
  
  
  require(rlang)
  require(magrittr)
  
  # Process GENE_RLOOP_OVERLAP
  message("Starting Genes")
  peakLst <- list(
    rloop_id = csvToGR(RLOOPS),
    gene_id = csvToGR(GENES)
  ) %>%
    makeOverlapTable() %>%
    readr::write_csv(file = GENE_RLOOP_OVERLAP)
  
  # Process GENFEAT_RLOOP_OVERLAP
  message("Starting GenFeat")
  peakLst <- list(
    rloop_id = csvToGR(RLOOPS),
    feature_id = csvToGR(GENFEATS)
  ) %>%
    makeOverlapTable() %>%
    readr::write_csv(file = GENFEAT_RLOOP_OVERLAP) 
  
  # xzip Compress
  message("Compress")
  system(paste0("xz -f ", GENE_RLOOP_OVERLAP))  
  system(paste0("xz -f ", GENFEAT_RLOOP_OVERLAP))
  
}

