#' Makes the "RLoopSignal" table for RMapDB
#'
#' @importFrom magrittr %>%
#' @import rlang
makeRLoopSignal <- function(RLOOPS, RLOOP_SIGNAL_MAT) {
  # Get RLoops
  rloops <- RLOOPS %>% csvToGR() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "names") %>%
    dplyr::mutate(loc = paste0(seqnames, "_", start, "_", end)) %>%
    dplyr::select(-c("seqnames", "start", "end", "width", "strand", "score"))
  
  # Read in the TSV
  readr::read_tsv(RLOOP_SIGNAL_MAT) %>% 
    dplyr::mutate(loc = paste0(`#'chr'`, "_", `'start'` + 1, "_", `'end'`)) %>%
    dplyr::select(-c(`#'chr'`, `'start'`,`'end'`)) %>%
    dplyr::right_join(rloops, by = "loc") %>%
    dplyr::select(-loc) %>%
    dplyr::rename(rloop_id=names) %>%
    tidyr::pivot_longer(cols = ! rloop_id,
                        names_to = "rmap_sample_id",
                        values_to = "norm_counts") %>%
    dplyr::mutate(rmap_sample_id = gsub(rmap_sample_id,
                                        pattern = "'([ES]+RX[0-9]+)_.+",
                                        replacement = "\\1"))
}


if (sys.nframe() == 0L | TRUE) {
  require(rlang)
  require(magrittr)
  source("analyses/Prepare-RMapDB-Tables/GeneRLoopOverlap_GenFeatRLoopOverlap/gene_rl_overlap__gf_rl_overlap.R")
  
  # TODO: Needs to become part of the snake pipe
  RLOOPS <- "analyses/Prepare-RMapDB-Tables/rloops.csv.xz"
  RMAPSAMPS <- "analyses/Prepare-RMapDB-Tables/rmap_samples.csv.xz"
  RLOOP_SIGNAL_MAT <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloop_signal.tsv.xz"
  RLOOP_SIGNAL_OUT <- "analyses/Prepare-RMapDB-Tables/rloop_signal.csv"
  
  # Need to create BED file for deeptools to use
  RLOOPS_BED <- "analyses/Prepare-RMapDB-Tables/RLoops/data/rloops.bed"
  RLOOPS %>% csvToGR() %>% rtracklayer::export.bed(con = RLOOPS_BED)
  system(paste0("xz -f ", RLOOPS_BED))
  
  # Make the RLoop Signal File
  message("Making RL signal output")
  makeRLoopSignal(RLOOPS = RLOOPS, RLOOP_SIGNAL_MAT = RLOOP_SIGNAL_MAT) %>%
    readr::write_csv(RLOOP_SIGNAL_OUT)
  
  # Compress
  message("Compress")
  system(paste0("xz -f ", RLOOP_SIGNAL_OUT))  
  
}

