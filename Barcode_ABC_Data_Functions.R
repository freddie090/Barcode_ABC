
# Copyright 2025 Cancer Research Technology and The Institute of Cancer Research.
#
# Licensed under a software academic use license provided with this software package (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at: https://github.com/freddie090/Barcode_ABC
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

################################################################################

# Description: a selection of data manipulation functions that process the 
# outputs of the Barcode_ABC simulations and inference.

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

#' Convert Wide Barcode Data Frame to Long Format
#'
#' This function converts a wide barcode data frame to long format while 
#' converting to relative frequencies and retaining only the top cumulative 
#' frequencies where >= crft.
#'
#' @param bc_df A wide-format data frame of barcode counts.
#' @param crft A numeric threshold for cumulative relative frequency.
#' @return A long-format data frame with relative frequencies and cumulative relative frequencies.
#' @export
bc_df_ltw <- function(bc_df, crft){
  # Convert counts to relative frequencies
  bc_df[, c(2:ncol(bc_df))] <- apply(bc_df[, c(2:ncol(bc_df))], 2, prop.table)
  # Turn from wide to long
  bc_dfl <- reshape2::melt(bc_df, id.vars = "bc", variable.name = "sample", value.name = "rf")
  # Split replicate and passage information into separate columns
  bc_dfl["rep"] <- sub("_P\\d{1}", "", bc_dfl$sample)
  bc_dfl["Passage"] <- sub("DT\\d{1}_", "", bc_dfl$sample)
  # Order by relative frequency and add a cumulative relative frequency column
  bc_dfl <- plyr::ddply(bc_dfl, .(rep, Passage), arrange, rf, decreasing = TRUE)
  bc_dfl <- plyr::ddply(bc_dfl, .(rep, Passage), function(x) {
    x$cum_rf <- cumsum(x$rf)
    x <- subset(x, cum_rf <= crft)
    return(x)
  })
  return(bc_dfl)
}

#' Custon Shannon Diversity function for the exception 'q=1' when calculating
#' # the Hill Numbers
#' @param xs A numeric vector of abundances
#' @param base 
shannon_diversity <- function(x, base = exp(1)) {
  # Ensure the input is numeric
  if (!is.numeric(x)) {
    stop("Input data must be numeric")
  }
  # Ensure the input data is non-negative
  if (any(x < 0, na.rm = TRUE)) {
    stop("Input data must be non-negative")
  }
  # Convert the input to a proportion
  x <- x / sum(x)
  # Calculate the Shannon diversity index
  shannon_index <- -sum(x * log(x, base = base), na.rm = TRUE)
  return(shannon_index)
}


#' Custom Diversity Index Function using Hill Numbers
#'
#' This function calculates a diversity index using Hill numbers.
#' It handles special cases for q = 1 and q = Inf.
#'
#' @param xs A numeric vector of abundances.
#' @param q A numeric value representing the diversity order.
#' @return The calculated diversity index.
#' @export
qD <- function(xs, q){
  xs <- xs[xs > 0]
  if(q == 1){
    #3 <- exp(diversity(xs, index = "shannon"))
    x3 <- exp(shannon_diversity(xs))
  } else if (q == Inf){
    x3 <- 1/(max(xs))
  } else {
    x1 <- xs^q
    x2 <- sum(x1)
    x3 <- x2^(1/(1 - q))
  }
  return(x3)
}


#' Calculate the qD Dissimilarity
#'
#' This function calculates the qD dissimilarity for a given data frame of 
#' abundances.
#'
#' @param xs_df A data frame where each column represents a sample and each row represents a lineage.
#' @param q A numeric value representing the diversity order.
#' @return The calculated qD dissimilarity.
#' @export
qD_diss <- function(xs_df, q){
  if(q == 1){
    # Collect log(xi/mean(x)), where
    # mean(x) = mean rel_freq x across all i's. 
    lxi_mx <- lapply(1:ncol(xs_df), function(i){
      log(xs_df[, i]/apply(xs_df, 1, mean))    
    })
    # Now want sum(xi * log(xi/mean(x)))
    lxi_mxxi <- sapply(1:ncol(xs_df), function(i){
      vals <- xs_df[, i]*lxi_mx[[i]]
      # Convert NAs to 0s
      vals[is.na(vals)] <- 0
      return(vals)
    })
    # Take sum across all x lineages
    mx1 <- apply(lxi_mxxi, 2, sum)
    # And now mean across all i populations
    mx2 <- mean(mx1)
    # And finally e^ this value  
    qM <- exp(mx2)
  } else if(q == Inf){
    # Here we only look at the most abundant
    # lineages. 
    max_xs <- xs_df[apply(xs_df, 2, which.max) ,]
    qM <- nrow(max_xs)
  } else {
    # Collect a qD | q for each sample. 
    qDs <- apply(xs_df, 2, qD, q=q)
    # Now all ^(1-q)
    qDs_pw <- qDs^(1 - q)
    # Now take the mean
    qDm <- mean(qDs_pw)
    # Now back-transform using ^1/(1-q)
    qDm_bt <- qDm^(1/(1 - q))
    # Calculate a qD for all pooled samples and
    # re-normalise so sum=1.0 for rel_freqs. 
    p_qDs <- rowSums(xs_df)/ncol(xs_df)
    # Now calculate pooled qD
    pqD <- qD(p_qDs, q=q)
    # Finally, want pooled/mean
    qM <- pqD/qDm_bt
  }
  return(qM)
}


#' Calculate Absolute Euclidean Distance in qD Statistic Space
#'
#' This function calculates the absolute Euclidean distance in qD statistic 
#' space.
#' The qD axis is transformed using log10(qD+1).
#'
#' @param qD1 First qD value.
#' @param qD2 Second qD value.
#' @param qD_diss1 First qD dissimilarity value.
#' @param qD_diss2 Second qD dissimilarity value.
#' @return The normalized Euclidean distance.
#' @export
qD_dist <- function(qD1, qD2, qD_diss1, qD_diss2,
                    Nmax=1e+08,
                    Nrep=4,
                    trans_qD=T,
                    trans_qD_diss=T){
                      
  if(trans_qD==T){
    qD1_trans <- (log10(qD1+1)/log10(Nmax+1))
    qD2_trans <- (log10(qD2+1)/log10(Nmax+1))
  }else{
    qD1_trans <- qD1
    qD2_trans <- qD2
  }

  if(trans_qD_diss==T){
    qD_diss1_trans <- ((qD_diss1 - 1)/(Nrep - 1))
    qD_diss2_trans <- ((qD_diss2 - 1)/(Nrep - 1))
  }else{
    qD_diss1_trans <- qD_diss1
    qD_diss2_trans <- qD_diss2
  }

  norm_qD_dist <- sqrt((qD1_trans - qD2_trans)^2 + (qD_diss1_trans - qD_diss2_trans)^2)
  return(norm_qD_dist)
}


#' Create a qD Statistic Data Frame
#'
#' This function takes a lineage data frame and returns a qD statistic data 
#' frame for a given passage.
#'
#' @param lin_df A data frame containing lineage data.
#' @param P A character string representing the passage.
#' @param q A numeric value representing the diversity order.
#' @return A data frame with columns for sample, qD, and qD dissimilarity.
#' @export
qD_df_fxn <- function(lin_df, P, q){
  temp_df <- lin_df[, grep(P, colnames(lin_df))]
  temp_df <- data.frame(apply(temp_df, 2, function(x){x/sum(x)}))
  qDs <- apply(temp_df, 2, qD, q=q)
  qD_disss <- qD_diss(temp_df, q=q)
  qD_df <- data.frame(samp = names(qDs), qD = qDs, qD_diss = qD_disss)
  return(qD_df)
}


#' Reorder Data Frame to Minimize qD Distance
#'
#' This function reorders a data frame of qD statistics to minimize the qD 
#' distance between all replicates given all possible iterations of simulated
#' vs ground-truth replicates.
#'
#' @param SIM_LIN_qD_DF Data frame of qD statistics for simulated lineages.
#' @param qD_GT_DF Data frame of ground truth qD statistics.
#' @param NREP Numeric value indiciating number of experimental replicates.
#' @param NPASS Numeric value indiciating number of experimental passages.
#' @param REP_NAs Logical indicating whether to replace NAs with maximum distance.
#' @param PLOT_VERS Logical indicating whether to return the plotting version.
#' @param DIST_VERS Numerical value indicating which qD fxn to use. Defaults to 1.
#' @param NMAX Numeric value indiciating how to normalise the qDs. Defaults to 1e+08.
#' @return A reordered data frame with minimized qD distance.
#' @export
min_qD_dist_df <- function(SIM_LIN_qD_DF, qD_GT_DF, NREP, NPASS,
                           REP_NAs = FALSE, PLOT_VERS = FALSE,
                           DIST_VERS = 1,
                           NMAX=1e+08){

  # Check that DIST_VERS is either 1 or 2
  if (!(DIST_VERS %in% c(1, 2))) {
    stop("DIST_VERS must be either 1 or 2.")
  }

  if(PLOT_VERS == FALSE){

      if(DIST_VERS==1){

        SAMP_PERMS <- data.frame(samp = paste0("DT", c(t(gtools::permutations(NREP, NREP, 1:NREP)))), 
                                gt_samp = paste0("DT", rep(1:NREP, factorial(NREP))),
                                set = rep(1:factorial(NREP), each = NREP))
        SAMP_PERMS <- bind_rows(replicate(NPASS, SAMP_PERMS, simplify = FALSE))
        SAMP_PERMS$P <- rep(as.character(1:NPASS), each = nrow(SAMP_PERMS)/NPASS)
        SIM_LIN_qD_DF$samp <- sub("_P\\d{1}", "", SIM_LIN_qD_DF$samp)
        qD_GT_DF$samp <- sub("_P\\d{1}", "", qD_GT_DF$samp)
        qD_GT_DF <- dplyr::rename(qD_GT_DF, gt_samp=samp)
        EXP_SIM_LIN_qD_DF <- left_join(SAMP_PERMS, SIM_LIN_qD_DF)
        EXP_SIM_LIN_qD_DF <- left_join(EXP_SIM_LIN_qD_DF, qD_GT_DF)

        EXP_SIM_LIN_qD_DF$qD_dist <- qD_dist(EXP_SIM_LIN_qD_DF$qD, 
                                            EXP_SIM_LIN_qD_DF$gt_qD, 
                                            EXP_SIM_LIN_qD_DF$qD_diss, 
                                            EXP_SIM_LIN_qD_DF$gt_qD_diss,
                                            Nmax=NMAX,
                                            Nrep=NREP)
      }
      if(DIST_VERS==2){
      # For use with the 'Re-Barcode' workflow - calculates distance solely on the 
      # new barcode 2nd Passage replicates.
        SAMP_PERMS <- data.frame(samp = paste0("DT", c(t(gtools::permutations(NREP, NREP, 1:NREP)))), 
                        gt_samp = paste0("DT", rep(1:NREP, factorial(NREP))),
                        set = rep(1:factorial(NREP), each = NREP))
        SAMP_PERMS$P <- "2"
        SIM_LIN_qD_DF$samp <- sub("_P\\d{1}", "", SIM_LIN_qD_DF$samp)
        qD_GT_DF$samp <- sub("_P\\d{1}", "", qD_GT_DF$samp)
        qD_GT_DF <- dplyr::rename(qD_GT_DF, gt_samp=samp)
        EXP_SIM_LIN_qD_DF <- left_join(SAMP_PERMS, SIM_LIN_qD_DF)
        EXP_SIM_LIN_qD_DF <- left_join(EXP_SIM_LIN_qD_DF, qD_GT_DF)

        EXP_SIM_LIN_qD_DF$qD_dist <- qD_dist(EXP_SIM_LIN_qD_DF$qD_2, 
                                             EXP_SIM_LIN_qD_DF$gt_qD_2,
                                             EXP_SIM_LIN_qD_DF$qD_diss_2,
                                             EXP_SIM_LIN_qD_DF$gt_qD_diss_2,
                                             Nmax=NMAX,
                                             Nrep=NREP)

      }

      if(REP_NAs == TRUE){

        EXP_SIM_LIN_qD_DF$qD_dist[is.na(EXP_SIM_LIN_qD_DF$qD_dist)] <- 1.0
        
      }

      SIM_LIN_SUMM <- plyr::ddply(EXP_SIM_LIN_qD_DF, .(set), summarise, 
                                  mean_qD_dist = mean(qD_dist),
                                  sd_qD_dist = sd(qD_dist))
                                  
      MIN_SET <- SIM_LIN_SUMM$set[first(which(SIM_LIN_SUMM$mean_qD_dist == min(SIM_LIN_SUMM$mean_qD_dist)))]
      FIN_LIN_qD_DF <- subset(EXP_SIM_LIN_qD_DF, set == MIN_SET)

  } else {

      EXP_SIM_LIN_qD_DF <- left_join(SIM_LIN_qD_DF, qD_GT_DF)

      if(DIST_VERS==1){
        EXP_SIM_LIN_qD_DF$qD_dist <- qD_dist(EXP_SIM_LIN_qD_DF$qD, 
                                            EXP_SIM_LIN_qD_DF$gt_qD, 
                                            EXP_SIM_LIN_qD_DF$qD_diss, 
                                            EXP_SIM_LIN_qD_DF$gt_qD_diss,
                                            Nmax=NMAX,
                                             Nrep=NREP)
      }
      if(DIST_VERS==2){
        EXP_SIM_LIN_qD_DF$qD_dist <- qD_dist(EXP_SIM_LIN_qD_DF$qD_2, 
                                             EXP_SIM_LIN_qD_DF$gt_qD_2,
                                             EXP_SIM_LIN_qD_DF$qD_diss_2,
                                             EXP_SIM_LIN_qD_DF$gt_qD_diss_2,
                                             Nmax=NMAX,
                                             Nrep=NREP)
      }

      FIN_LIN_qD_DF <- EXP_SIM_LIN_qD_DF
  }
  return(FIN_LIN_qD_DF)
}


#' Fill Missing Replicates in Data Frame
#'
#' This function inserts artificial counts for every replicate that doesn't 
#' have a barcode count.
#'
#' @param SIM_LIN_DF A data frame of simulated lineages.
#' @return A data frame with missing replicates filled.
#' @export
fill_empty_reps <- function(SIM_LIN_DF, NREP, NPASS){
  REP_NAMES <- paste0("DT", 1:NREP)
  REP_P_NAMES <- as.vector(outer(REP_NAMES, 
                                 paste0("_P", 1:NPASS), 
                                 FUN = paste0))
  MISS_REPS <- setdiff(REP_P_NAMES, colnames(SIM_LIN_DF))
  OUT_SIM_LIN_DF <- SIM_LIN_DF
  if(length(MISS_REPS)>0){
    for(i in seq_along(MISS_REPS)){
      OUT_SIM_LIN_DF[[MISS_REPS[[i]]]] <- 0
    }
  }
  return(OUT_SIM_LIN_DF)
}


#' Generate Data Frame of Posterior Distribution
#'
#' This function returns a data frame of the posterior distribution for a 
#' given parameter.
#'
#' @param POST_DF A data frame of posterior samples.
#' @param SYMBOL_DF A data frame of parameter symbols.
#' @param CHOSEN_PAR A character string representing the chosen parameter.
#' @param USE_GT A logical value indicating whether to include ground truth parameters. Defaults to FALSE.
#' @param GT_PARAMS A data frame with the ground truth parametes. Empty by default.
#' @return A data frame of posterior samples for the chosen parameter.
#' @export
prod_post_df <- function(POST_DF, SYMBOL_DF, CHOSEN_PAR,
                         USE_GT=F, GT_PARAMS=NA){
  PAR_POST_DF <- POST_DF[, c(CHOSEN_PAR, "sim_id")]
  PAR_POST_DF <- gather(PAR_POST_DF, key = "param", value = "val", 1)
  PAR_POST_DF$param <- as.factor(PAR_POST_DF$param)
  if(USE_GT==T){
    PAR_POST_DF["gt_val"] <- subset(GT_PARAMS, param == CHOSEN_PAR)$gt_val
  }
  levels(PAR_POST_DF$param) <- SYMBOL_DF$symbols[match(levels(PAR_POST_DF$param), SYMBOL_DF$param)]
  return(PAR_POST_DF)
}


#' Read CSV files and add samp_n column
#'
#' @param file A character string representing the file path.
#' @return A data frame with an added column samp_n extracted from the file name.
read_and_add_samp_n <- function(file) {
  df <- read.csv(file)
  samp_n <- as.numeric(sub(".*samp_(\\d*).csv", "\\1", file))
  df <- dplyr::mutate(df, samp_n = samp_n)
  return(df)
}

#' Read CSV files and add samp_n column - 2nd version. 
#'
#' @param file A character string representing the file path.
#' @return A data frame with an added column samp_n extracted from the file name.
read_and_add_samp_n_2 <- function(file) {
  df <- read.csv(file)
  samp_n <- as.numeric(sub(".*_df_(\\d*).csv", "\\1", file))
  df <- dplyr::mutate(df, samp_n = samp_n)
  return(df)
}


#' Read and combine CSV files for a specific generation
#'
#' @param gen A character string representing the generation identifier.
#' @param prefix A character string representing the prefix for the file names.
#' @return A combined data frame of all the CSV files matching the generation and prefix.
read_generation_files <- function(gen, prefix) {
  pattern <- paste0(prefix, "_gen_", gen, "_post")
  files <- grep(pattern, list.files(), value = TRUE)
  dfs <- purrr::map(files, read_and_add_samp_n)
  dplyr::bind_rows(dfs)
}


#' Read and combine CSV files for a specific generation
#' Version for use with the Plasticity Modulator simulation workflow.
#' 
#' @param gen A character string representing the generation identifier.
#' @param prefix A character string representing the prefix for the file names.
#' @param suffix A character string representing the suffix for the filenames.
#' @return A combined data frame of all the CSV files matching the generation and prefix.
read_generation_files_plast_mod <- function(gen, prefix, suffix) {
  pattern <- paste0(prefix, "_gen_", gen, "_post_", suffix)
  files <- grep(pattern, list.files(), value = TRUE)
  dfs <- purrr::map(files, read_and_add_samp_n)
  dplyr::bind_rows(dfs)
}


#' Read in the PPC files for cell size trajectory plot.
#' 
#' @param sim_code A charactter - 1 or 2 - denoting which simulation type. 
#' @param df_type A character string taking the value 'pop' or 'sol' indicating the type of dataframe being read in. 
#' @return A combined data frame of all the CSV files matching the dataframe type.
read_PPC_pop_or_sol_dfs <- function(sim_code, df_type){
  if(!(sim_code %in% c("1", "2"))){
    stop("sim_code should be '1' or '2'.")
  }
  if(!(df_type %in% c("pop", "sol"))){
    stop("df_type should be 'pop' or 'sol'.")
  }
  pattern <- paste0("PPC_", sim_code, "_ode_", df_type, "_df\\d+.csv")
  files <- grep(pattern, list.files(), value = TRUE)
  dfs <- purrr::map(files, read_and_add_samp_n_2)
  df <- dplyr::bind_rows(dfs)
  if(df_type == "pop"){
    df <- df %>%
      dplyr::select(t,u,samp_n) %>%
      dplyr::rename("N" = "u")
  }
  return(df)
}


#' Load posterior distributions for a specific generation
#'
#' @param gen A character string representing the generation identifier.
#' @param turn_long A boolean determining whether to convert to a long dataframe.
#' Default = FALSE.
#' @return A data frame of the posterior distribution with columns renamed and gathered.
load_posterior_distribution <- function(file_regexp, turn_long = FALSE) {
  file <- grep(file_regexp, list.files(), value = TRUE)
  post_df <- read.csv(file) %>%
    dplyr::rename("del" = "delt") %>%
    dplyr::mutate(del = 10^-(del))
  if(turn_long==T){
    post_df <- post_df  %>%
      tidyr::gather(key = "param", value = "val", 2:ncol(post_df)) %>%
      dplyr::select(param, val)
  }
  return(post_df)
}


#' Process ground truth parameters based on model type
#'
#' @param param_file A character string representing the path to the parameters CSV file.
#' @param mod_type A character string representing the simulation type (e.g., "A").
#' @param param_vec A vector of parameter names to keep.
#' @param symbol_df A data frame with parameter symbols.
#' @return A data frame with processed ground truth parameters.
#' @export
process_ground_truth_params <- function(param_file, mod_type, 
                                        param_vec, symbol_df){

  gt_params <- read.csv(param_file)

  gt_params$rho <- -log10(gt_params$rho)
  gt_params$mu <- -log10(gt_params$mu)
  
  if (mod_type %in% c("B", "C")) {
    gt_params$sig <- -log10(gt_params$sig)
  }
  if (mod_type == "C") {
    gt_params$al <- -log10(gt_params$al)
  }
  
  gt_params <- gt_params[, param_vec]
  gt_params <- tidyr::gather(gt_params, key = "param", value = "gt_val")
  
  gt_params$sym_param <- as.factor(gt_params$param)
  levels(gt_params$sym_param) <- symbol_df$symbols[match(levels(gt_params$sym_param), symbol_df$param)]
  
  return(gt_params)

}


#' Process ground truth parameters based on model type
#' This version is specific to the 'Plasticity Modulator' simulation workflow.
#'
#' @param param_file A character string representing the path to the parameters CSV file.
#' @param mod_type A character string representing the model type (e.g., "A").
#' @param param_vec A vector of parameter names to keep.
#' @param symbol_df A data frame with parameter symbols.
#' @param sim_code A character string representing which Plast.Mod simulation: 1 or 2.
#' @return A data frame with processed ground truth parameters.
#' @export
process_ground_truth_params_plast_mod <- function(param_file, mod_type, 
                                                  param_vec, symbol_df,
                                                  sim_code){

  gt_params <- read.csv(param_file)

  # Sim 1 - assumes constant single mu throughout.
  if(sim_code=="1"){
    gt_params <- gt_params %>%
      dplyr::rename("mu" = "mu_1") %>%
      dplyr::select(-c(mu_2))
      gt_params$mu <- -log10(gt_params$mu)
  }
  # Sim 2 - different mu for pre- (mu_1) and during- treatment (mu_2).
  if(sim_code=="2"){
    gt_params$mu_1 <- -log10(gt_params$mu_1)
    gt_params$mu_2 <- -log10(gt_params$mu_2)
  }

  gt_params$rho <- -log10(gt_params$rho)
  
  if (mod_type %in% c("B", "C")) {
    gt_params$sig <- -log10(gt_params$sig)
  }
  if (mod_type == "C") {
    gt_params$al <- -log10(gt_params$al)
  }
  
  gt_params <- gt_params[, param_vec]
  gt_params <- tidyr::gather(gt_params, key = "param", value = "gt_val")
  
  gt_params$sym_param <- as.factor(gt_params$param)
  levels(gt_params$sym_param) <- symbol_df$symbols[match(levels(gt_params$sym_param), symbol_df$param)]
  
  return(gt_params)

}


#' Load and process observed data
#'
#' @param sol_file A character string representing the path to the solution data CSV file.
#' @param pop_file A character string representing the path to the population data CSV file.
#' @return A list containing two data frames: processed solution data and processed population data.
#' @export
load_and_process_observed_data <- function(sol_file, pop_file, nrep) {
  # Load the data
  gt_sol_df <- read.csv(sol_file)
  gt_pop_df <- read.csv(pop_file)
  
  # Add replicate information and subset the data
  gt_pop_df["rep"] <- rep(unique(gt_sol_df$rep), each = nrow(gt_pop_df) / nrep)
  gt_pop_df <- subset(gt_pop_df, t != 4.0)
  gt_pop_df <- gt_pop_df[, c("t", "u", "rep")]
  colnames(gt_pop_df) <- c("t", "N", "rep")
  
  return(list(sol_df = gt_sol_df, pop_df = gt_pop_df))
}


#' Process ground truth lineage data
#'
#' @param lin_file A character string representing the path to the lineage data CSV file.
#' @param q A numeric value representing the diversity order. Defaults to q=2.
#' @return A data frame with processed ground truth qD values for P1 and P2 populations.
#' @export
process_ground_truth_lin_df <- function(lin_file, q=2) {
  # Read the ground truth lineage data
  gt_lin_df <- read.csv(lin_file)
  
  # Calculate qD values for P1 and P2 populations
  qD_gt_P1_df <- qD_df_fxn(gt_lin_df, P = "P1", q = q)
  qD_gt_P1_df["P"] <- "1"
  qD_gt_P2_df <- qD_df_fxn(gt_lin_df, P = "P2", q = q)
  qD_gt_P2_df["P"] <- "2"
  
  # Combine the results and rename columns
  qD_gt_df <- bind_rows(qD_gt_P1_df, qD_gt_P2_df)
  colnames(qD_gt_df) <- sub("qD", "gt_qD", colnames(qD_gt_df))
  
  return(qD_gt_df)
}


#' Read and preprocess wide barcode dataframe
#'
#' @param file_path A character string representing the file path to the CSV file.
#' @return A preprocessed data frame with modified column names.
#' @export
read_and_preprocess_barcode_data <- function(file_path) {
  df <- read.csv(file_path)
  colnames(df)[1] <- "bc"
  colnames(df) <- sub("_ce1_cluster.csv", "", colnames(df))
  return(df)
}


#' Postprocess diversity metrics dataframe
#'
#' @param df A data frame containing diversity metrics.
#' @return A postprocessed data frame with updated column names.
#' @export
process_gt_diversity_dataframe <- function(df, P, q, remove_regexp) {
  result_df <- qD_df_fxn(df, P = P, q = q)
  result_df$P <- as.character(P)
  result_df$P <- sub("P", "", result_df$P)
  colnames(result_df) <- sub("qD", "gt_qD", colnames(result_df))
  result_df$samp <- sub(remove_regexp, "", result_df$samp)
  return(result_df)
}


#' Read CSV files matching a pattern
#'
#' @param pattern A character string representing the pattern to match files.
#' @return A list of data frames read from the matching CSV files.
#' @export
read_csv_files_by_pattern <- function(pattern) {
  files <- grep(pattern, list.files(), value = TRUE)
  dfs <- sapply(files, read.csv, simplify = FALSE, USE.NAMES = TRUE)
  return(dfs)
}


#' Extract qD statistics for a given population
#'
#' @param dfs A list of data frames containing simulation data.
#' @param P A character string representing the population identifier.
#' @param q A numeric value representing the diversity order.
#' @return A list of data frames with qD statistics.
#' @export
extract_sim_qD_statistics <- function(dfs, P, q) {
  qD_dfs <- lapply(dfs, qD_df_fxn, P = P, q = q)
  return(qD_dfs)
}


#' Add simulation IDs to data frames
#'
#' @param dfs A list of data frames.
#' @param pattern A character string representing the pattern to extract sim IDs.
#' @return A list of data frames with added sim_id column.
#' @export
add_simulation_ids <- function(dfs, pattern) {
  for(i in seq_along(dfs)) {
    dfs[[i]]["sim_id"] <- sub(pattern, "\\1", names(dfs)[[i]])
  }
  return(dfs)
}


#' Postprocess diversity statistics data frame
#'
#' @param df A data frame containing diversity statistics.
#' @return A postprocessed data frame with updated row names and added population column.
#' @export
postprocess_diversity_stats <- function(df) {
  rownames(df) <- 1:nrow(df)
  df["P"] <- sub(".*_P(\\d{1})", "\\1", df$samp)
  return(df)
}


#' Load and preprocess parameter dataframe
#'
#' @param file_path A character string representing the file path to the CSV file.
#' @param param_vec A vector of parameter names to keep.
#' @return A preprocessed data frame with the selected parameters.
#' @export
load_and_preprocess_param_df <- function(file_path, param_vec) {
  df <- read.csv(file_path)
  if("delt" %in% colnames(df)){
  df <- df %>%
    dplyr::rename("del" = "delt")
  } 
  df <- df %>%
    dplyr::select(param_vec)
  df["sim_id"] <- as.character(1:nrow(df))
  return(df)
}

#' Subset and transform parameter data frame for plotting
#'
#' @param param_df A data frame containing the parameter dataframe
#' @param lin_stat_df A data frame containing the diversity statistics.
#' @param symbol_df A data frame with parameter symbols.
#' @param subset_by_sim_id A logical value indicating whether to subset by sim_id. Defaults to TRUE.
#' @return A transformed data frame suitable for plotting.
#' @export
transform_param_df_for_plotting <- function(param_df, lin_stat_df, symbol_df, subset_by_sim_id = TRUE) {
  if (subset_by_sim_id) {
    df <- param_df %>%
      subset(sim_id %in% unique(lin_stat_df$sim_id)) %>%
      dplyr::select(-c(sim_id)) %>%
      gather(key = "param", value = "val")
  } else {
    df <- param_df %>%
      dplyr::select(-c(sim_id)) %>%
      gather(key = "param", value = "val")
  }
  df$param <- as.factor(df$param)
  levels(df$param) <- symbol_df$symbols[match(levels(df$param), symbol_df$param)]
  return(df)
}


#' Calculate observed vector of summary statistics
#'
#' @param gt_qD_df A data frame containing ground truth qD values.
#' @return A numeric vector of summary statistics.
#' @export
get_obsv_qD_vec <- function(gt_qD_df) {
  gt_qD_df$rep <- sub("DT(\\d+)_P\\d+", replacement = "\\1", gt_qD_df$samp)
  gt_qD_df$P <- sub("DT\\d+_P(\\d+)", replacement = "\\1", gt_qD_df$samp)
  gt_qD_df <- dplyr::arrange(gt_qD_df, P, rep)
  gt_sstats <- c(gt_qD_df$gt_qD, gt_qD_df$gt_qD_diss)
  return(gt_sstats)
}


#' Calculate observed vector of summary statistics
#' A second version that works with the re-labelled qD statistics in the 
#' 'Re-Barcode' simulation workflow.
#'
#' @param gt_qD_df A data frame containing ground truth qD values.
#' @return A numeric vector of summary statistics.
#' @export
get_obsv_qD_vec_2 <- function(gt_qD_df) {
  gt_qD_df$rep <- sub("DT(\\d+)_P\\d+", replacement = "\\1", gt_qD_df$samp)
  gt_qD_df$P <- sub("DT\\d+_P(\\d+)", replacement = "\\1", gt_qD_df$samp)
  gt_qD_df <- dplyr::arrange(gt_qD_df, P, rep)
  gt_sstats <- c(gt_qD_df$gt_qD_2, gt_qD_df$gt_qD_diss_2)
  return(gt_sstats)
}


#' Returns a dataframe of the qD statistics ordered by passage and replicate.
#'
#' @param lin_stat_df A data frame containing linear statistics.
#' @param nrep A numberical value denoting the number of replicates.
#' @return A data frame with NAs replaced and replicates labeled.
#' @export
get_sim_qD_vec <- function(lin_stat_df, nrep) {
  lin_stat_df <- plyr::ddply(lin_stat_df, .(sim_id, P), function(x) {
    x$rep <- rep(1:nrep)
    return(x)
  })
  
  lin_sim_stats <- plyr::ddply(lin_stat_df, .(sim_id), function(x) {
    dplyr::select(x, c(rep, P, qD, qD_diss)) %>%
      tidyr::pivot_wider(names_from = P, values_from = c(qD, qD_diss)) %>%
      tidyr::pivot_wider(names_from = rep, values_from = dplyr::contains("qD"))
  }) %>%
    dplyr::select(-c(sim_id))
  
  return(lin_sim_stats)
}


#' Returns a dataframe of the qD statistics ordered by passage and replicate.
#' A second version that works with the re-labelled qD statistics in the 
#' 'Re-Barcode' simulation workflow.
#'
#' @param lin_stat_df A data frame containing linear statistics.
#' @param nrep A numberical value denoting the number of replicates.
#' @return A data frame with NAs replaced and replicates labeled.
#' @export
get_sim_qD_vec_2 <- function(lin_stat_df, nrep) {
  lin_stat_df <- plyr::ddply(lin_stat_df, .(sim_id, P), function(x) {
    x$rep <- rep(1:nrep)
    return(x)
  })
  
  lin_sim_stats <- plyr::ddply(lin_stat_df, .(sim_id), function(x) {
    dplyr::select(x, c(rep, P, qD_2, qD_diss_2)) %>%
      tidyr::pivot_wider(names_from = P, values_from = c(qD_2, qD_diss_2)) %>%
      tidyr::pivot_wider(names_from = rep, values_from = dplyr::contains("qD"))
  }) %>%
    dplyr::select(-c(sim_id))
  
  return(lin_sim_stats)
}


#' Normalize qD Values
#'
#' Normalizes `qD` values using a log10 transformation based on a maximum value.
#'
#' @param qD A numeric vector of `qD` values to normalize.
#' @param NMAX A numeric value representing the maximum normalisation factor.
#' @return A numeric vector of normalized `qD` values.
#' @examples
#' norm_qD(c(10, 100, 1000), NMAX = 1e8)
norm_qD <- function(qD, NMAX) {
  nqD <- log10(qD + 1) / log10(NMAX + 1)
  return(nqD)
}

# Vectorize the `norm_qD` function for element-wise operations
norm_qD <- Vectorize(norm_qD)

#' Normalize qD Dissimilarity Values
#'
#' Normalizes `qD_diss` values based on the number of replicates.
#'
#' @param qD_diss A numeric vector of dissimilarity `qD` values to normalize.
#' @param NREP An integer specifying the number of replicates.
#' @return A numeric vector of normalized dissimilarity `qD` values.
#' @examples
#' norm_qD_diss(c(1.5, 2.0, 3.0), NREP = 4)
norm_qD_diss <- function(qD_diss, NREP) {
  nqD_diss <- (qD_diss - 1) / (NREP - 1)
  return(nqD_diss)
}

# Vectorize the `norm_qD_diss` function for element-wise operations
norm_qD_diss <- Vectorize(norm_qD_diss)

#' Calculate Deviance
#'
#' Computes the deviance between observed summary statistics and posterior samples.
#'
#' @param observed_ss A numeric vector of observed summary statistics.
#' @param posterior_ss A matrix or data frame of posterior summary statistics (rows correspond to samples).
#' @return A numeric value representing the deviance.
#' @examples
#' observed_ss <- c(0.2, 0.4, 0.6)
#' posterior_ss <- matrix(runif(30), nrow = 10)
#' deviance_fun(observed_ss, posterior_ss)
deviance_fun <- function(observed_ss, posterior_ss) {
  ss_dist <- apply(posterior_ss, 1, function(row) sqrt(sum((row - observed_ss)^2)))
  ss_epsi <- max(ss_dist)
  ss_dist <- ss_dist / ss_epsi
  ss_dist <- exp(-0.5 * (ss_dist^2)) / sqrt(2 * pi)
  ss_dist <- log(ss_dist / ss_epsi)
  deviance <- -2 * mean(ss_dist)
  return(deviance)
}

#' Calculate Deviance Information Criterion (DIC)
#'
#' Computes the DIC for model evaluation using observed and posterior summary statistics.
#'
#' @param observed_ss_path A string specifying the file path to the observed summary statistics CSV.
#' @param post_ss_path A string specifying the file path to the posterior summary statistics CSV.
#' @param post_m_ss_path A string specifying the file path to the posterior mean summary statistics CSV.
#' @param NMAX A numeric value representing the maximum normalisation factor for `qD`.
#' @param NREP An integer specifying the number of replicates for `qD_diss`.
#' @return A numeric value representing the DIC.
#' @examples
#' # Assuming `observed.csv`, `posterior.csv`, and `posterior_mean.csv` exist:
#' calculate_DIC("observed.csv", "posterior.csv", "posterior_mean.csv", NMAX = 1e8, NREP = 4)
calculate_DIC <- function(observed_ss_path,
                          post_ss_path,
                          post_m_ss_path,
                          NMAX,
                          NREP) {
  # Read and normalize observed data
  observed_qD <- read.csv(observed_ss_path) %>%
    dplyr::pull(gt_qD) %>%
    norm_qD(NMAX = NMAX)
  
  observed_qD_diss <- read.csv(observed_ss_path) %>%
    dplyr::pull(gt_qD_diss) %>%
    norm_qD_diss(NREP = NREP)
  
  observed_ss <- c(observed_qD, observed_qD_diss)
  
  # Read and normalize posterior data
  ppc_sstats <- read.csv(post_ss_path) %>%
    dplyr::mutate(across(-contains("_diss_"), ~ norm_qD(.x, NMAX = NMAX))) %>%
    dplyr::mutate(across(contains("_diss_"), ~ norm_qD_diss(.x, NREP = NREP)))
  
  pD_sstats <- read.csv(post_m_ss_path) %>%
    dplyr::mutate(across(-contains("_diss_"), ~ norm_qD(.x, NMAX = NMAX))) %>%
    dplyr::mutate(across(contains("_diss_"), ~ norm_qD_diss(.x, NREP = NREP)))
  
  # Compute DIC
  D_bar <- deviance_fun(observed_ss, ppc_sstats)
  D.theta_bar <- deviance_fun(observed_ss, pD_sstats)
  DIC <- 2 * D_bar - D.theta_bar
  
  return(DIC)
}


################################################################################