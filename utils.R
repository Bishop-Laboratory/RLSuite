# Taken from https://github.com/cran/VennDiagram/blob/master/R/hypergeometric.test.R
calculate.overlap.and.pvalue = function(list1, list2, total.size, lower.tail = TRUE, adjust = FALSE) {
  
  # calculate actual overlap
  actual.overlap <- length(intersect(list1, list2));
  
  # calculate expected overlap
  # need to cast to avoid integer overflow when length(list1) * length(list2) is extremely large
  expected.overlap <- as.numeric(length(list1)) * length(list2) / total.size;
  
  adjust.value <- 0;
  
  # adjust actual.overlap to reflect P[X >= x]
  if (adjust & !lower.tail) {
    adjust.value <- 1;
    warning('Calculating P[X >= x]');
  }
  
  # calculate significance of the overlap
  overlap.pvalue <- phyper(
    q = actual.overlap - adjust.value,
    m = length(list1),
    n = total.size - length(list1),
    k = length(list2),
    lower.tail = lower.tail
  );
  
  # return values
  return( c(actual.overlap, expected.overlap, overlap.pvalue) );
  
}

plot_single_feature <- function(feature, db, pltdat, lmts, yb) {
  
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
plot_feature <- function(feature, db, pltdat, lmts, yb, splitby) {
  
  # Define limits
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
        x = splitby,
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
    ggplot2::xlab(paste0("RLSeq ", toupper(splitby))) +
    ggplot2::labs(
      title = paste0(feature, " - ", db),
      subtitle = "Feature enrichment results"
    ) +
    ggpubr::stat_compare_means(comparisons = list(unique(pltdat[,splitby, drop=TRUE])),
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

plot_multi_feature <- function(features, db, pltdat, lmts, yb, factororder, axislmts) {
  DODGE_POS <- .8
  data4plot <-  pltdat %>%
    as.list() %>%
    bind_rows() %>%
    filter(type %in% {{ features }}, 
           db %in% {{ db }}) %>%
    filter(stat_fisher_rl > lmts[1], stat_fisher_rl < lmts[2]-1,
           ! is.na(prediction)) %>%
    mutate(type = factor(type, levels = factororder),
           cond = paste0(label, "_", prediction),
           cond = factor(cond, levels = c("NEG_NEG", "POS_NEG", "NEG_POS", "POS_POS")))
  kw <- data4plot %>% 
    rstatix::group_by(type) %>%
    rstatix::kruskal_test(stat_fisher_rl ~ cond) 
  print(kw)
  
  df_p_val <- data4plot %>% 
    rstatix::group_by(type) %>%
    rstatix::dunn_test(stat_fisher_rl ~ cond, p.adjust.method = "bonferroni") %>%
    rstatix::add_significance(p.col = "p") %>% 
    rstatix::add_xy_position(dodge = DODGE_POS, x = "type")
  df_p_val <- df_p_val %>%
    filter(group1 == "NEG_NEG")
  
  # Define limits
  splitby <- "prediction"
  cols <- RLSeq:::auxdata[[paste0(tolower(splitby[1]), "_cols")]]
  cols <- tribble(
    ~cond, ~col,
    "NEG_NEG", "#8a2c2c",
    "POS_NEG", "#A76767",
    "NEG_POS", "#7aa4c4",
    "POS_POS", "#2270ab"
  )
  colvec <- cols %>% dplyr::pull(.data$col)
  names(colvec) <- as.data.frame(cols)[, "cond"]
  yval <- "stat_fisher_rl"
  data4plot %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = "type",
        y = yval,
        label = "experiment",
        pattern_alpha="cond",
        pattern_fill="label",
        fill = "cond",
        color = "cond"
      )
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_violin(
      color = "black",
      width = 1,
      trim = FALSE,
      alpha = .5,
      na.rm = TRUE,
      position = ggplot2::position_dodge(.8)
    ) +
    geom_boxplot(
      width = .3,
      color = "black",
      alpha = 1,
      position = ggplot2::position_dodge(.8)
    ) +
    ggplot2::scale_fill_manual(
      values = setNames(cols$col, nm = cols$cond)
    ) +
    ggplot2::ylab(
      "Fisher's test odds ratio (log2)"
    ) +
    ggplot2::xlab("RLSeq Prediction") +
    ggplot2::labs(
      title = db,
      subtitle = "Feature enrichment results"
    ) +
    theme_prism(base_size = 14) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 14),
      plot.caption = ggplot2::element_text(size = 11)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(size = 0, stroke = 0)
      )
    ) +
    guides(fill = guide_legend(title="Condition")) +
    scale_y_continuous(limits = axislmts) +
    add_pvalue(df_p_val,
               xmin = "xmin",
               xmax = "xmax",
               inherit.aes = F)
}

plot_multi_feature_iptype <- function(features, db, pltdat, lmts, factororder) {
  DODGE_POS <- .8
  data4plot <-  pltdat %>%
    as.list() %>%
    bind_rows() %>%
    filter(type %in% {{ features }}, 
           db %in% {{ db }}) %>%
    filter(stat_fisher_rl > lmts[1], stat_fisher_rl < lmts[2]-1,
           ! is.na(ip_type)) %>%
    mutate(type = factor(type, levels = factororder),
           cond = factor(ip_type))
  df_p_val <- data4plot %>% 
    rstatix::group_by(type) %>%
    rstatix::wilcox_test(stat_fisher_rl ~ cond)  %>%
    rstatix::add_significance(p.col = "p") %>% 
    rstatix::add_xy_position(dodge = DODGE_POS, x = "type")
  print(df_p_val)
  
  # Define limits
  cols <- tribble(
    ~cond, ~col,
    "dRNH", "#bf79b2",
    "S9.6", "#30a7d6"
  )
  colvec <- cols %>% dplyr::pull(.data$col)
  names(colvec) <- as.data.frame(cols)[, "cond"]
  yval <- "stat_fisher_rl"
  data4plot %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        x = "type",
        y = yval,
        label = "experiment",
        pattern_alpha="cond",
        pattern_fill="cond",
        fill = "cond",
        color = "cond"
      )
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_violin(
      color = "black",
      width = 1,
      trim = FALSE,
      alpha = .5,
      na.rm = TRUE,
      position = ggplot2::position_dodge(.8)
    ) +
    geom_boxplot(
      width = .3,
      color = "black",
      alpha = 1,
      position = ggplot2::position_dodge(.8)
    ) +
    ggplot2::scale_fill_manual(
      values = setNames(cols$col, nm = cols$cond)
    ) +
    ggplot2::ylab(
      "Fisher's test odds ratio (log2)"
    ) +
    ggplot2::xlab("Histone Marks") +
    ggplot2::labs(
      title = db,
      subtitle = "Feature enrichment results"
    ) +
    theme_prism(base_size = 14) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 14),
      plot.caption = ggplot2::element_text(size = 11)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(size = 0, stroke = 0)
      )
    ) +
    ggpubr::rotate_x_text(45, vjust = 1) +
    guides(fill = guide_legend(title="Condition")) +
    # scale_y_continuous(limits = axislmts) +
    add_pvalue(df_p_val,
               xmin = "xmin",
               xmax = "xmax",
               inherit.aes = F)
}


# Deeptools metaplot
dtools_metaplot <- function(db, type) {
  matfile <- paste0("data/annot_mats/", db,"__", type, ".mat.gz")
  metap <- read_tsv(matfile, skip = 1, col_names = FALSE)
  metadat <- read_lines(matfile, n_max = 1)
  metadat <- gsub(metadat, pattern = "@", replacement = "") %>% jsonlite::fromJSON()
  mat <- metap %>%
    select(-c(1:6)) 
  mat %>%
    colMeans(na.rm = TRUE) -> cm
  divs <- setNames(c(42, 184), nm = c("rlregions_dRNH", "rlregions_S96"))
  pltdats <- tibble("score" = cm,
                    "loc" = rep(seq(600), 2),
                    "div" = c(rep(divs[metadat$sample_labels[1]], metadat$sample_boundaries[2]),
                              rep(divs[metadat$sample_labels[2]], metadat$sample_boundaries[2])),
                    "label" = c(rep(metadat$sample_labels[1], metadat$sample_boundaries[2]),
                                rep(metadat$sample_labels[2], metadat$sample_boundaries[2]))) %>%
    group_by(label) %>%
    mutate(score = score / div,
           #score = scale(score, center = F),
           # score = scales::rescale(x = score, to = c(0, 1)),
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
      breaks = c(1, 300, 600),
      labels = c("-3000", type, "+3000"),
      expand = c(0,0)
    ) +
    ggpubr::rremove("legend") +
    scale_y_continuous(expand = c(0,0),
                       limits = c(mcm, mxcm),
                       breaks = mcm + c(0*mcdiff, .25*mcdiff, .5*mcdiff, .75*mcdiff, 1*mcdiff),
                       labels = c(0, .25, .5, .75, 1)) +
    geom_area(aes(group = label), alpha = .5, position = "identity")  +
    ggtitle(paste0(db, " - ", type))
    return(plt)
}





