
# Barcode_ABC_Plotting_Functions.R
# Description: a selection of plotting functions that plot the outputs of the 
# Barcode_ABC simulations and inference.

# Load required libraries 
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(purrr)

################################################################################

# Function Definitions


#' Custom Plotting Theme
#'
#' This function creates a custom ggplot2 theme for barcode plots.
#'
#' @return A ggplot2 theme object.
#' @param text_size Size of text for Barcode plotting theme. Defaults to 16. 
#' @export
theme_BARCODE <- function(text_size=16){
  theme_minimal() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = text_size))
}

#' Plot Color Swatch
#'
#' This function plots a color swatch for a given vector of colors.
#'
#' @param x A vector of colors (hex, numeric, or string).
#' @return A color swatch plot.
#' @export
swatch <- function(x) {
  par(mai=c(0.2, max(strwidth(x, "inch") + 0.4, na.rm = TRUE), 0.2, 0.4))
  barplot(rep(1, length(x)), col=rev(x), space = 0.1, axes=FALSE,
          names.arg=rev(x), cex.names=0.8, horiz=TRUE, las=1)
}

#' Generate Colors Using iwanthue Algorithm
#'
#' This function generates a specified number of colors using the iwanthue 
#' algorithm.
#'
#' @param n Number of colors.
#' @param hmin Lower bound of hue (0-360).
#' @param hmax Upper bound of hue (0-360).
#' @param cmin Lower bound of chroma (0-180).
#' @param cmax Upper bound of chroma (0-180).
#' @param lmin Lower bound of luminance (0-100).
#' @param lmax Upper bound of luminance (0-100).
#' @param plot Should a color swatch be plotted? (default: FALSE).
#' @param random Should clustering be random? (default: FALSE).
#' @return A vector of generated colors.
#' @export
iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100,
                     plot=FALSE, random=FALSE) {
  require(colorspace)
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0,
            hmax <= 360, cmax <= 180, lmax <= 100,
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            n > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1),
                                   seq(-100, 100, 5),
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin &
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    lab <- as(hcl, 'LAB')
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), n, iter.max=50)
  if (isTRUE(plot)) {
    swatch(hex(LAB(clus$centers)))
  }
  hex(LAB(clus$centers))
}


#' Custom Label Maker for Plotting
#'
#' This function creates custom labels for plotting. It formats specific Greek letters (e.g., \U03C1, \U03BC, \U03C3, \U03B1, \U03B4) using `label_bquote`.
#'
#' @param value A character vector of values to be labelled.
#' @return A character vector of formatted labels.
#' @export
make_label <- function(value) {
  x <- as.character(value)
  if (x %in% c("\U03C1", "\U03BC", "\U03C3", "\U03B1", "\U03B4")) {
    x <- label_bquote(10 ^ - (.x))
    return(x)
  } else {
    return(x)
  }
}


#' Custom Plot Labeller
#'
#' This function creates custom labels for ggplot2 facets using the `make_label` function.
#'
#' @param variable The name of the variable to be labelled.
#' @param value The value of the variable to be labelled.
#' @return A character vector of formatted labels suitable for `ggplot2`.
#' @export
plot_labeller <- function(variable, value) {
  do.call(expression, lapply(levels(value), make_label))
}


#' Plot Observed vs Posterior Total Population Size Change Through Treatment
#'
#' This function plots the observed vs posterior total population size change 
#' through treatment.
#'
#' @param TREAT_DF A data frame of treatment information.
#' @param SOL_DF A data frame of solution information.
#' @param POP_DF A data frame of population information.
#' @param POST_LABEL A character string representing the posterior label.
#' @return A ggplot2 object of the population trajectories.
#' @export
plot_obsv_post_pop_traj <- function(TREAT_DF, SOL_DF, POP_DF, POST_LABEL, 
                                    NMAX, TMAX){
  p <- ggplot() +
    geom_rect(data = TREAT_DF, 
              aes(ymin = 1.0, ymax = NMAX, 
                  xmin = start, xmax = end,
                  fill = treat), alpha = 0.2, colour = "black") +
    scale_fill_manual(values = c("palegreen1", "red"), name = "Treatment") + 
    geom_line(data = SOL_DF, aes(x = t, y = N, 
                                 group = interaction(rep, samp_n),
                                 colour = type),
              alpha = 0.2, size = 0.8) +
    geom_point(data = subset(POP_DF, type == POST_LABEL), 
               aes(x = t, y = N, colour =  POST_LABEL), size = 3) + 
    geom_point(data = subset(POP_DF, type == "Ground Truth"), 
               aes(x = t, y = N, colour =  "Ground Truth"), size = 3) + 
    scale_colour_manual(values = c("black", "blue"), name = POST_LABEL) + 
    scale_y_log10(limits = c(1, NMAX)) + 
    scale_x_continuous(limits = c(0, TMAX)) + 
    xlab("Time (days)") +
    ylab("N\n") +
    theme_BARCODE() +
    ggtitle(paste0("Population Trajectories: ", POST_LABEL))
  return(p)
}


#' Plot Lineage Statistics as Function of Parameter
#'
#' This function plots the lineage statistics (qD) as a function of one of the 
#' given parameters.
#'
#' @param LIN_STAT_QD_DF_PV A data frame of lineage statistics for plotting.
#' @param LIN_STAT_QD_DF A data frame of lineage statistics.
#' @param QD_DIST_SUMM_DF A data frame of qD distance summaries.
#' @param CHOSEN_PAR A character string representing the chosen parameter.
#' @param PLOT_GT_VAL A logical value indicating whether to plot the ground-truth parameter value. Defaults to FALSE.
#' @param GT_PARAMS A DataFrame with the ground-truth parameters. Defaults to empty.
#' @return A ggplot2 object of the lineage statistics plot.
#' @export
lin_stat_plot <- function(LIN_STAT_QD_DF_PV, LIN_STAT_QD_DF, 
                          QD_DIST_SUMM_DF,
                          CHOSEN_PAR,
                          PLOT_GT_VAL=F,
                          GT_PARAMS=NA,
                          PAR_DF=NA,
                          SAMP_PAR_DF=NA,
                          TEXT_SIZE=16){
  
  PAR_LIM_L <- min(subset(param_lims, key == CHOSEN_PAR)$val)
  PAR_LIM_U <- max(subset(param_lims, key == CHOSEN_PAR)$val)
  PAR_LIM_RNG <- (PAR_LIM_U - PAR_LIM_L)
  PAR_SYM <- unique(subset(param_lims, key == CHOSEN_PAR)$param)
  
  p1 <- ggplot(data = subset(LIN_STAT_QD_DF_PV, qD_type == "qD"), 
               aes(x = qD_val+1, y = qD_diss_val)) + 
    geom_line(data = LIN_STAT_QD_DF_PV, 
              aes(x = qD_val+1, y = qD_diss_val, 
                  group = interaction(samp, sim_id))) + 
    geom_point(aes_string(fill = CHOSEN_PAR),
               shape = 21, colour = "black", size = 5, alpha = 0.8) + 
    geom_point(data = subset(LIN_STAT_QD_DF_PV, qD_type == "gt_qD"),
               aes(x = qD_val+1, y = qD_diss_val), colour = "red", size = 4) +
    scale_x_log10(limits = c(1e+00, 1e+06),
                  breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + 
    scale_y_continuous(limits = c(1, 4)) + 
    xlab(expression(""^{"q=2"}~D)) + 
    ylab(expression(""^{"q=2"}~D~({"\U03B2"}))) +
    scale_fill_viridis_c(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) +
    facet_grid(P~rep) + 
    theme_BARCODE(text_size=TEXT_SIZE) +
    theme(panel.grid.major = element_line(size = 0.5, colour = "grey90"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.0),
          legend.position = "none")
  
  p2 <- ggplot(data = QD_DIST_SUMM_DF, 
               aes_string(x = CHOSEN_PAR, y = "qD_dist_mean", 
                          fill = CHOSEN_PAR)) +
    geom_jitter(shape = 21, size = 5, alpha = 0.6,
                width = PAR_LIM_RNG/30) +
    scale_fill_viridis_c(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) +
    scale_colour_viridis_c(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) +
    scale_x_continuous(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) + 
    scale_y_continuous(limits = c(0, 0.8)) +
    ylab(expression(atop(""^{"q=2"}~D, paste("Distance")))) +
    theme_BARCODE(text_size=TEXT_SIZE)
  
  if(PLOT_GT_VAL==T){
    p2 <- p2 + 
      geom_vline(xintercept = subset(GT_PARAMS, param == CHOSEN_PAR)$gt_val,
                 linetype = "dashed", size = 1, colour = "red")
  }
  
  LIN_STAT_PLOT <- cowplot::plot_grid(p1, p2, nrow = 1, 
                                      rel_widths = c(1.0, 0.5), 
                                      rel_heights = c(1.0, 0.7))
  return(LIN_STAT_PLOT)
}


#' Plot Lineage Statistics as Function of Parameter
#'
#' This function plots the lineage statistics (qD) as a function of one of the 
#' given parameters - this 2nd version only uses the within-replicate diversity
#' # (qD) for use with the Re-Barcoded simulation outputs. 
#'
#' @param LIN_STAT_QD_DF_PV A data frame of lineage statistics for plotting.
#' @param LIN_STAT_QD_DF A data frame of lineage statistics.
#' @param QD_DIST_SUMM_DF A data frame of qD distance summaries.
#' @param CHOSEN_PAR A character string representing the chosen parameter.
#' @param PLOT_GT_VAL A logical value indicating whether to plot the ground-truth parameter value. Defaults to FALSE.
#' @param GT_PARAMS A DataFrame with the ground-truth parameters. Defaults to empty.
#' @return A ggplot2 object of the lineage statistics plot.
#' @export
lin_stat_plot_2 <- function(LIN_STAT_QD_DF_PV, LIN_STAT_QD_DF, 
                            QD_DIST_SUMM_DF,
                            CHOSEN_PAR,
                            PLOT_GT_VAL=F,
                            GT_PARAMS=NA,
                            TEXT_SIZE=16){
  
  PAR_LIM_L <- min(subset(param_lims, key == CHOSEN_PAR)$val)
  PAR_LIM_U <- max(subset(param_lims, key == CHOSEN_PAR)$val)
  PAR_LIM_RNG <- (PAR_LIM_U - PAR_LIM_L)
  PAR_SYM <- unique(subset(param_lims, key == CHOSEN_PAR)$param)
  
  p1 <- ggplot(data = subset(LIN_STAT_QD_DF_PV, qD_type == "qD_2"), 
              aes_string(y = "qD_2_val+1", x = CHOSEN_PAR)) + 
    # geom_line(data = LIN_STAT_QD_DF_PV, 
    #           aes_string(y = "qD_2_val+1", x = CHOSEN_PAR, 
    #               group = "interaction(samp, sim_id)")) + 
    geom_point(aes_string(fill = CHOSEN_PAR),
              shape = 21, colour = "black", size = 5, alpha = 0.8) + 
    geom_point(data = subset(LIN_STAT_QD_DF_PV, qD_type == "gt_qD_2"),
              aes_string(x = "qD_2_val+1", y = CHOSEN_PAR), size = 4) +
    scale_y_log10(limits = c(1e+00, 1e+05),
                  breaks = scales::trans_breaks("log10", function(x) 10^x), 
                  labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + 
    scale_x_continuous(limits = c(PAR_LIM_L, PAR_LIM_U)) + 
    ylab(expression(""^{"q=2"}~D)) + 
    scale_fill_viridis_c(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) +
    facet_grid(P~rep) + 
    theme_BARCODE(text_size=TEXT_SIZE) +
    theme(panel.grid.major = element_line(size = 0.5, colour = "grey90"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.0),
          legend.position = "none")
  
  p2 <- ggplot(data = QD_DIST_SUMM_DF, 
               aes_string(x = CHOSEN_PAR, y = "qD_dist_mean", 
                          fill = CHOSEN_PAR)) +
    geom_jitter(shape = 21, size = 5, alpha = 0.6,
                width = PAR_LIM_RNG/30) +
    # geom_jitter(data = LIN_STAT_QD_DF, aes_string(x = CHOSEN_PAR, y = "qD_dist"),
    #             shape = 21, size = 4, alpha = 0.6,
    #             width = PAR_LIM_RNG/30) +
    scale_fill_viridis_c(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) +
    scale_colour_viridis_c(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) +
    scale_x_continuous(limits = c(PAR_LIM_L, PAR_LIM_U), name = PAR_SYM) + 
    scale_y_continuous(limits = c(0, 0.3)) +
    ylab(expression(atop(""^{"q=2"}~D, paste("Distance")))) +
    theme_BARCODE(text_size=TEXT_SIZE)
  
  if(PLOT_GT_VAL==T){

    GT_PARAMS_qD_2_DF <- LIN_STAT_QD_DF_PV %>% 
      dplyr::filter(qD_type=="gt_qD_2") %>% 
      dplyr::select(qD_2_val, rep) %>% 
      unique() %>% 
      dplyr::mutate(x=subset(gt_params, param==CHOSEN_PAR)$gt_val)

    p1 <- p1 + 
    geom_point(data = GT_PARAMS_qD_2_DF,
               aes_string(y = "qD_2_val+1", x = "x"),
               colour = "red", size = 4) #+
      # geom_hline(yintercept = subset(GT_PARAMS, param == CHOSEN_PAR)$gt_val,
      #            linetype = "dashed", size = 1, colour = "red") 
    p2 <- p2 + 
      geom_vline(xintercept = subset(GT_PARAMS, param == CHOSEN_PAR)$gt_val,
                 linetype = "dashed", size = 1, colour = "red")
  }
  
  LIN_STAT_PLOT <- cowplot::plot_grid(p1, p2, nrow = 1, 
                                      rel_widths = c(1.0, 0.5), 
                                      rel_heights = c(1.0, 0.7))
  return(LIN_STAT_PLOT)
}


#' Plot and save posterior distribution as boxplots for a single parameter
#'
#' @param chosen_par A character string representing the chosen parameter.
#' @param abc_df A data frame containing the ABC posterior data.
#' @param param_lims A data frame containing the parameter limits.
#' @param plot_gt_val A logical value indicating whether to plot the ground-truth parameter value. Defaults to FALSE.
#' @return A ggplot object representing the posterior distribution boxplot for the chosen parameter.
#' @export
plot_post_box <- function(chosen_par, abc_df, param_lims,
                          plot_gt_val=FALSE,
                          col_vec=col_vec) {
  
  generate_labels <- function(powers) {
    sapply(powers, function(x) as.expression(bquote(10^.(x))))
  }
  
  post_box <- ggplot(data = subset(abc_df, param == chosen_par), 
                     aes(x = chosen_par, y = val)) + 
    geom_blank(data = param_lims[param_lims$key %in% param_vec ,], 
               aes(y = val)) + 
    geom_boxplot(fill = col_vec[4], alpha = 0.5) + 
    scale_y_continuous(limits = c(
      min(subset(param_lims, param == chosen_par)$val),
      max(subset(param_lims, param == chosen_par)$val))) +
    theme_BARCODE() + 
    xlab("") +
    ylab("")
  
  if(plot_gt_val==T){
    post_box <- post_box + 
      geom_hline(aes(yintercept = gt_val), size = 1, 
                 linetype = "dashed", colour = "red")
  }
  
  if(chosen_par %in% c("mu", "mu_1", "mu_2", "sig", "al")) {
    post_box <- post_box + 
      scale_y_reverse(limits = c(
        max(subset(param_lims, param == chosen_par)$val), 
        min(subset(param_lims, param == chosen_par)$val)),
        breaks = c(9, 7, 5, 3, 1),
        labels = generate_labels(c(-9, -7, -5, -3, -1)),
        name = "")
  }
  
  if(chosen_par %in% c("rho")) {
    post_box <- post_box + 
      scale_y_reverse(limits = c(
        max(subset(param_lims, param == chosen_par)$val), 
        min(subset(param_lims, param == chosen_par)$val)),
        breaks = c(0, 2, 4, 6, 8),
        labels = generate_labels(c(0, -2, -4, -6, -8)),
        name = "")
  }
  
  return(post_box)
}

################################################################################