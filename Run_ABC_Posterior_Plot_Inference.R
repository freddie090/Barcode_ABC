
# Takes the full ABC-lineage inference output files from a single 
# parameter/observation set, calculates the qD diversity statistics and 
# creates summary output plots illustrating the inference process and 
# calculates the final posterior distribution by retaining the simulations
# with the closest diversity statistics.

rm(list = ls())

################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) == 9 || length(args) == 10)) {
  stop("There should be 9 or 10 command-line arguments.")
}

# Load in command line arguments:
abc_scr_dir <- args[1]
abc_mod_dir <- args[2]
abc_out_dir <- args[3]
abc_model <- args[4]
abc_fit_dir <- args[5]
output_type <- args[6]
is_simulated_data <- args[7]
param_df_name <- args[8]
calculate_posterior <- as.logical(args[9])
if (calculate_posterior && length(args) == 10) {
  post_n <- as.integer(args[10])
}

################################################################################

# Experimental Parameters:
##########################

# Treatment Timings

treat_ons <- seq(4.0, 100.0, 8.0)
treat_ons <- rep(treat_ons, each = 2)
treat_offs <- seq(8.0, 100.0, 8.0)
treat_offs <- rep(treat_offs, each = 2)
treat_offs <- c(0, treat_offs)

length(treat_ons) <- min(length(treat_ons), length(treat_offs))
length(treat_offs) <- min(length(treat_ons), length(treat_offs))

treat_df <- data.frame("start" = treat_offs, "end" = treat_ons)
treat_df["treat"] <- rep(c("off", "on"), length.out = nrow(treat_df))

# Number of replicates:
nrep <- 4

# Number of Passages:
npass <- 2

################################################################################
 
# Load in custom plotting and data analysis functions:
source(paste0(abc_scr_dir, "/Barcode_ABC_Data_Functions.R"))
source(paste0(abc_scr_dir, "/Barcode_ABC_Plotting_Functions.R"))

################################################################################

# Load the parameter limits into a dataframe to set the facet_wrap limits

param_lims <- data.frame(
  param = c(rep("Dc", 2), 
            rep("del", 2), 
            rep("k", 2), 
            rep("psi", 2), 
            rep("rho", 2),
            rep("mu", 2),
            rep("sig", 2),
            rep("al", 2)),
  val = c(c(0, 8),
          c(0, 1),
          c(0, 4),
          c(0, 1),
          c(0, 7), 
          c(0, 9), 
          c(0, 9),
          c(0, 9))
)

# Parameter vector given the simulation type: 

if(abc_model == "A"){
  param_vec <- c("Dc", "del", "k", "psi", "rho", "mu")
}
if(abc_model == "B"){
  param_vec <- c("Dc", "del", "k", "psi", "rho", "mu", "sig")
}
if(abc_model == "C"){
  param_vec <- c("Dc", "del", "k", "psi", "rho", "mu", "sig", "al")
}

# Create a dataframe of symbols used in plotting with corresponding character
# names. 

symbol_df <- data.frame(param = c("rho", "mu", "sig", "al", "del", 
                                  "psi", "Dc", "k"),
                        symbols = c("\U03C1", "\U03BC", "\U03C3", "\U03B1", "\U03B4", 
                                    "\U03C8", "Dc", "\U03BA"))

sym_param_vec <- symbol_df[symbol_df$param %in% param_vec ,]$symbols

# Convert character parameters to symbols for plotting, keeping the original 
# character parameter names as keys: 

param_lims$key <- param_lims$param
param_lims$param <- as.factor(param_lims$param)
levels(param_lims$param) <- symbol_df$symbols[match(levels(param_lims$param), symbol_df$param)]

# Create a colour vector for plotting: 

col_vec <- colorRampPalette(brewer.pal("Paired", n = 12))(12)

################################################################################

# Modify param lims to set facet limits via 'geom_blank': 

param_lims_mod <- bind_rows(param_lims, param_lims)
param_lims_mod["P"] <- rep(c(1, 2), each = nrow(param_lims))
param_lims_mod["sim_id"] <- "0"
param_lims_mod["qD_dist"] <- 0
param_lims_mod["qD_dist_mean"] <- 0
param_lims_mod["param_dist"] <- 0
param_lims_mod["keep"] <- T


# Set tmax and Nmax for population size change plots:

plot_tmax <- 80.0
plot_Nmax <- 1e+08

################################################################################

# GO to output directory:

setwd(abc_out_dir)
setwd(abc_fit_dir)


# Load in the Ground Truth Parameters if Simulated Data:
########################################################

if(is_simulated_data == T){
  
  param_file <- "param_df.csv"
  
  gt_params <- process_ground_truth_params(param_file, abc_model,
                                           param_vec, symbol_df)
  
}


# Ground Truth Lineage Distribution
###################################

# Define the path to the lineage data CSV file
lin_file <- "lin_lin_df.csv"

# Process the ground truth lineage data
qD_gt_df <- process_ground_truth_lin_df(lin_file)


# Simulated Lineage Distributions
#################################

# Read all of the simulated barcode dataframes from the agent-based 
# simulations: 

sim_lin_dfs <- read_csv_files_by_pattern(paste0("KMC_", output_type, 
                                                "_lin_df_[0-9].*.csv"))

# Calculate the qD statistics for the simulated barcode distributions:  

sim_qD_df <- lapply(paste0("P", 1:npass),
                    extract_sim_qD_statistics,
                    dfs=sim_lin_dfs,
                    q=2) %>%
  unlist(recursive = F) %>%
  add_simulation_ids(pattern = ".*_df_(\\d.*).csv") %>%
  bind_rows() %>%
  postprocess_diversity_stats()

# Load and preprocess the simulation parameters dataframe:

param_df <- load_and_preprocess_param_df(param_df_name, param_vec)

# Transform posterior data frame for plotting, subsetting by sim_id
param_dfl <- transform_param_df_for_plotting(param_df, sim_qD_df, symbol_df)

# Transform posterior data frame for plotting, including entire posterior distribution
param_full_dfl <- transform_param_df_for_plotting(param_df, sim_qD_df,
                                                  symbol_df,
                                                  subset_by_sim_id = FALSE)

# Add to the simulated ABC dataframe
sim_qD_df <- left_join(sim_qD_df, param_df, by = "sim_id")

# Transform del to [0, 1]: 
sim_qD_df$del <- 10^-(sim_qD_df$del)

# Calculate distance in qD Statistic space:
###########################################

# Do a plotting version that doesn't re-order the replicates
# to minimise the distance: 
sim_qD_df_pv <- plyr::ddply(sim_qD_df, .(sim_id), 
                            min_qD_dist_df, qD_GT_DF=qD_gt_df,
                            NREP=nrep, NPASS=npass,
                            PLOT_VERS=T)

# And a version that re-orders to minimise qD statistics for finding the 
# closest simulations:
sim_qD_df <- plyr::ddply(sim_qD_df, .(sim_id), 
                         min_qD_dist_df, qD_GT_DF=qD_gt_df,
                         NREP=nrep, NPASS=npass,
                         REP_NAs = F, PLOT_VERS=F)

# Turn to long via parameter column 
sim_qD_dfl <- gather(sim_qD_df, key = "param", value = "val", 
                     grep(paste(param_vec, collapse = "|"), 
                          colnames(sim_qD_df)))

# Join the ground-truth param dataframe onto lineage statistic dataframe
sim_qD_dfl <- left_join(sim_qD_dfl, gt_params, by = c("param"))

# Set parameter symbols
sim_qD_dfl$param <- as.factor(sim_qD_dfl$param)
levels(sim_qD_dfl$param) <- symbol_df$symbols[match(levels(sim_qD_dfl$param), symbol_df$param)]


# Process qD dataframes for plotting

lin_stat_qD <- sim_qD_df_pv %>% 
  pivot_longer(ends_with("qD"), names_to = "qD_type", 
               values_to = "qD_val") %>%
  dplyr::select(-c("qD_diss", "gt_qD_diss"))

lin_stat_qDdiss <- sim_qD_df_pv %>%
  pivot_longer(ends_with("qD_diss"), names_to = "qD_type",
               values_to = "qD_diss_val") %>%
  dplyr::select(-c("qD", "gt_qD"))

lin_stat_qDdiss$qD_type <- sub("_diss", "", lin_stat_qDdiss$qD_type)

sim_qD_dfl_pv <- join(lin_stat_qD, lin_stat_qDdiss)

sim_qD_dfl_pv["rep"] <- sub("_P\\d{1}", "", sim_qD_dfl_pv$samp)
sim_qD_dfl_pv["P"] <- paste0("P", sim_qD_dfl_pv$P)

# Can change the number of simulations to be plotted:
plot_nsim <- 20

sim_qD_dfl_pv_sub <- filter(sim_qD_dfl_pv,
                            sim_id %in% (sim_qD_dfl_pv$sim_id %>%
                                           unique() %>%
                                           sample(plot_nsim, replace = F)))

sim_qD_df_sub <- sim_qD_df[sample(1:nrow(sim_qD_df),
                                  plot_nsim*nrep,
                                  replace = F) ,]

###############################################################################

######################################
# Plotting ABC Population Trajectories
######################################


# Simulated population and solution dataframes:
KMC_pop_df <- grep(paste0("KMC_", output_type, "_pop_df"), list.files(), value = T) %>%
    purrr::map(read_and_add_samp_n_2) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(type = paste0("ABC_", output_type))

KMC_sol_df <- grep(paste0("KMC_", output_type, "_sol_df"), list.files(), value = T) %>%
    purrr::map(read_and_add_samp_n_2) %>%
    dplyr::bind_rows() %>%
  dplyr::mutate(type = paste0("ABC_", output_type))


# Ground truth dataframes:
pop_df <- read.csv("lin_pop_df.csv") %>%
  dplyr::mutate(type = "Ground Truth")
KMC_pop_df <- bind_rows(KMC_pop_df, pop_df) %>%
  dplyr::rename(N = u)

if(is_simulated_data==T){
  sol_df <- read.csv("lin_sol_df.csv") %>%
    dplyr::mutate(type = "Ground Truth")
  KMC_sol_df <- bind_rows(KMC_sol_df, sol_df)
}


ABC_pop_traj_plot <- plot_obsv_post_pop_traj(treat_df, 
                                             KMC_sol_df, KMC_pop_df, 
                                             paste0("ABC_", 
                                                    output_type),
                                             plot_Nmax, plot_tmax)

################################################################################

# Calculate Posterior
#####################

# Create a summary dataframe of the mean diversity distance per parameter
# set (simulation id): 
sim_qD_summ_df <- plyr::ddply(sim_qD_df, .(sim_id), 
                              summarise, 
                              qD_dist_mean = mean(qD_dist)) %>%
  dplyr::left_join(sim_qD_df %>%
                     dplyr::select(any_of(c(param_vec, "sim_id"))) %>%
                     unique()) %>%
  dplyr::arrange(qD_dist_mean)


if(calculate_posterior==T){
  
  # Just keep the (post_tol * nsim) closest simulations according to their 
  # mean qD distance - opposed to the 'abc' package normalising the values
  # first:
  
  abc_df <- sim_qD_summ_df[1:post_n, param_vec]
  
  # Save the posterior dataframe
  
  write.csv(abc_df, 
            file = "ABC_post_dist.csv",
            row.names = F)
  
  abc_df <- abc_df %>% 
    gather(key = "param", value = "val")
  
  abc_df <- left_join(abc_df, gt_params)
  abc_df$param <- abc_df$sym_param 
  
  plot_ncol <- switch(abc_model,
                      "A" = 3,
                      "B" = 4,
                      "C" = 4)
  
  # Plot final posterior boxplot:
  abc_post_box_plot <- lapply(sym_param_vec, plot_post_box, 
                              abc_df = abc_df,
                              param_lims = param_lims,
                              plot_gt_val = is_simulated_data,
                              col_vec=col_vec)
  
  abc_post_box_plot <- cowplot::plot_grid(plotlist=abc_post_box_plot,
                                          ncol = plot_ncol)
  
}

################################################################################

# Save Plots
############

if(dir.exists("plots")==F){
  dir.create("plots")
  setwd("plots")
} else {
  setwd("plots")
}

# Change ggsave function so bg = white
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)



# Posterior Plot
################

if(calculate_posterior==T){
  
  plot_width <- switch(abc_model,
                       "A" = 8,
                       "B" = 10,
                       "C" = 10)
  
  ggsave("ABC_posterior_plot.jpg",
         plot = abc_post_box_plot,
         height = 11,
         width = plot_width)
  
}


# Cell Population Change Plot
#############################

ggsave(paste0("ABC_", output_type, "_pop_traj_plot.jpg"), 
       ABC_pop_traj_plot, 
       height = 4, 
       width = 12)


# qD Statistic Plots
####################

# Function to create plots for parameters
create_plots <- function(sim_qD_dfl_pv_sub, sim_qD_df_sub, 
                         qD_dist_summ_df, symbols,
                         gt_params) {
  plots <- mapply(function(sym) {
    lin_stat_plot(LIN_STAT_QD_DF_PV=sim_qD_dfl_pv_sub, 
                  LIN_STAT_QD_DF=sim_qD_df_sub, 
                  QD_DIST_SUMM_DF=qD_dist_summ_df,
                  CHOSEN_PAR=sym,
                  PLOT_GT_VAL=T,
                  GT_PARAMS=gt_params)
  }, symbols, SIMPLIFY = FALSE)
  names(plots) <- symbols
  plots
}


# Function to save plots
save_plots <- function(plots, symbols) {
  lapply(symbols, function(sym) {
    ggsave2(paste0("ABC_", output_type, "_lineage_qD_", sym, "_plot.jpg"),
            plot = plots[[sym]],
            height = 6,
            width = 17)
  })
}

# Define symbols for each simulation type
symbols_list <- list(
  A = c("rho", "mu", "del"),
  B = c("rho", "mu", "del", "sig"),
  C = c("rho", "mu", "del", "sig", "al")
)

# Create parameter and sample parameter data frames based on simulation type
if(abc_model %in% names(symbols_list)) {
  symbols <- symbols_list[[abc_model]]
  
  # Create plots
  plots <- create_plots(sim_qD_dfl_pv_sub, sim_qD_df_sub, 
                        sim_qD_summ_df, symbols, 
                        gt_params)
  
  # Save plots
  save_plots(plots, symbols)
}

################################################################################

