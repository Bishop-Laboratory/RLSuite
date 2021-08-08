#' Makes the "Genomic Features" table for RMapDB
#'
#' @importFrom magrittr %>%
#' @import rlang
makeGenomicFeatures <- function() {
    annotatr::build_annotations(genome = "hg38", annotations = c(
      "hg38_cpg_islands", "hg38_basicgenes", "hg38_lncrna_gencode", 
      "hg38_enhancers_fantom", "hg38_genes_promoters", "hg38_genes_intergenic"
    )) %>%
    as.data.frame() %>%
    dplyr::mutate(
      location = paste0(
        seqnames, ":", start, "-", end, ":", strand
      )
    ) %>%
    dplyr::select(
      id, type, source=tx_id, location
    )
}

if (interactive()) {
  
  GENOMIC_FEATURES <- "analyses/Prepare-RMapDB-Tables/genomic_features.csv"
  
  require(rlang)
  require(magrittr)
  
  readr::write_csv(makeGenomicFeatures(), file = GENOMIC_FEATURES)
  system(paste0("xz ", GENOMIC_FEATURES))  # xzip
  
}
