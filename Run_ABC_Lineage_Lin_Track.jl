
using RCall
using DataFrames

# Simulates synthetic lineage sizes through treatment given known parameter 
# values. Uses the agent-based kinetic monte carlo (KMC) model to simulate 
# individual barcode lineages, retaining the top x lineages at periodic times
# throughout the simulations. These lineages are plotted alongside the diversity
# statistics at each Passage time point at the end of the simulation (as in 
# Figure 1 in the manuscript). 
# Simulation type (Model A - Preex, Model B - Switch, Model C - Escape) and 
# 'ground truth' parameter values are loaded from respective command line 
# arguments.

length(ARGS) == 5 || error("Should be 5 command line arguments.")

abc_mod_dir = ARGS[1]
abc_out_dir = ARGS[2]
abc_model = ARGS[3]
param_table_path = ARGS[4]
sim_id = Int64(eval(Meta.parse(ARGS[5])))

abc_model ∈ ["A", "B", "C"] || error("The 'abc_model' variable is not one of the three permitted values: A, B or C.")

# Load in both sets of functions for the phenotypic compartment ('Population')
# and agent-based kinetic monte-carlo (KMC) models. 
include(string(abc_mod_dir, "/ABC_Population_Model_Fxns.jl"))
include(string(abc_mod_dir, "/ABC_Barcode_KMC_Model.jl"))

# Read in the parameter table
ptbl_df = DataFrame(CSV.File(string(param_table_path)))

# Print current time to output
println(now())
# Check number of threads being used.
println(Base.Threads.nthreads())

cd(abc_out_dir)

# Load in the variable 'ground truth' parameters from the command line. 

if abc_model == "A"

    rho=Float64(ptbl_df[sim_id, :rho])
    mu=Float64(ptbl_df[sim_id, :mu])
    Dc=Float64(ptbl_df[sim_id, :Dc])
    del=Float64(ptbl_df[sim_id, :del])
    k=Float64(ptbl_df[sim_id, :k])
    psi=Float64(ptbl_df[sim_id, :psi])

    sig = 0.0
    al = 0.0

end

if abc_model == "B"

    rho=Float64(ptbl_df[sim_id, :rho])
    mu=Float64(ptbl_df[sim_id, :mu])
    sig=Float64(ptbl_df[sim_id, :sig])
    Dc=Float64(ptbl_df[sim_id, :Dc])
    del=Float64(ptbl_df[sim_id, :del])
    k=Float64(ptbl_df[sim_id, :k])
    psi=Float64(ptbl_df[sim_id, :psi])

    al = 0.0

end

if abc_model == "C"

    rho=Float64(ptbl_df[sim_id, :rho])
    mu=Float64(ptbl_df[sim_id, :mu])
    sig=Float64(ptbl_df[sim_id, :sig])
    al=Float64(ptbl_df[sim_id, :al])
    Dc=Float64(ptbl_df[sim_id, :Dc])
    del=Float64(ptbl_df[sim_id, :del])
    k=Float64(ptbl_df[sim_id, :k])
    psi=Float64(ptbl_df[sim_id, :psi])

end

################################################################################

########################################################
# Set the fixed parameters for the simulated experiment:
######################################################## 

# Number of uniquely barcoded cells when the experiment begins:
n0 = Int64(1e+04);

# Estimates of the average, sensitive population's birth and death rates:
b = 0.893; d = 0.200;

# The measurement error of the population size readings (where
# observed N ~ Normal(mean = true N, sd = true N * phi)):
phi = 0.10;

# Time when the experiment begins and maximum time of the experiment (in days):
t0 = 0.0; tmax = 60.0;

# Maximum number of cells before a replicate is harvested, and an estimate of
# the carrying capacity:
Nmax = Int64(80*1e+04); Cc = Int64(100*1e+04);

# The treatment windows when treatment begins (treat_ons) and ends (treat_offs)
# (in days):
treat_ons = collect(4.0:8.0:tmax);
treat_offs = collect(8.0:8.0:tmax);

# The time to Passage cells (in days) and number of Passages (currently 
# max n_Pass = 2):
t_Pass = 30.0; n_Pass = 2

# The expansion time of the barcoded cells before splitting into experimental
# replicates (in days):
t_exp = 6.0;

# The observation times to record the population size changes (in days):
t_keep = [5.0, 10.0];

# The time-window when calculating the drug-concentration change, and a Boolean
# whether to treat the cells (leave these unless troubleshooting):
dt_save_at = 1e-03; treat = true; 

# The population size per-phenotype to switch from a stochastic jump process
# (when <= Nswitch) to a deterministic ODE approximation (> Nswitch).
Nswitch = 100;

# The number of highest count barcoded lineages to retain throughout the 
# simulation (recorded at tmax/t_frac).
top_x = 10;

################################################################################

###############
# Run the Model
###############

lin_sim_out = exp_pulse_treat_kmc_full_sim_fxn(n0, b, d, rho, mu, sig, del, al,
                                                Dc, k, psi,
                                                t0, t_exp, tmax, t_Pass,
                                                Nmax, Cc,
                                                treat_ons, treat_offs,
                                                dt_save_at, t_keep,
                                                n_Pass=n_Pass, R_real="l",
                                                top_x=top_x, lin_track=true)

# Barcode lineage observed populations and parameters
lin_pop_df = DataFrame(t = lin_sim_out["t"], u = lin_sim_out["u"])

# Barcode lineage full solution population trajectories
lin_sol_df = lin_sim_out["sol_df"]

# Barcode lineage distribution 
lin_lin_df = lin_sim_out["lin_df"]

# Top x lineages over time
sub_lin_df = vcat(lin_sim_out["sub_lin_dfs"]...)

# This simulation's parameters 
param_df = DataFrame(abc_model = abc_model, 
                     rho = rho, mu = mu, sig = sig, al = al, del = del,                        
                     Dc = Dc, k = k, psi = psi, phi = phi,
                     n0 = n0, b = b, d = d, Nmax = Nmax, Cc = Cc)

# Use RCall to save dataframes and generate plots.
@rput lin_pop_df; @rput lin_lin_df; 
@rput lin_sol_df;
@rput sub_lin_df; 
@rput param_df;
@rput treat_ons; @rput treat_offs; @rput tmax;
@rput abc_model; @rput sim_id;

R"""

    # Assumes helper functions are stored in parental directory: 
    source("../Barcode_ABC_Data_Functions.R")
    source("../Barcode_ABC_Plotting_Functions.R")

    # For sub-sampling sampled times when plotting: 

    find_closest <- function(value, vector) {
    abs_diff <- abs(vector - value)
    return(which.min(abs_diff))
    }

    # Make Model and simulation directories
    if(dir.exists(paste0("Model_", abc_model, "_outputs"))==F){
        dir.create(paste0("Model_", abc_model, "_outputs"))
        setwd(paste0("Model_", abc_model, "_outputs"))
    } else {
        setwd(paste0("Model_", abc_model, "_outputs"))
    }
    if(dir.exists(paste0("sim_", sim_id))==F){
        dir.create(paste0("sim_", sim_id))
        setwd(paste0("sim_", sim_id))
    } else {
        setwd(paste0("sim_", sim_id))
    }

    # Save dataframes
    write.csv(lin_pop_df, "lin_pop_df.csv", row.names = F)
    write.csv(lin_lin_df, "lin_lin_df.csv", row.names = F)
    write.csv(lin_sol_df, "lin_sol_df.csv", row.names = F)
    write.csv(sub_lin_df, "sub_lin_df.csv", row.names = F)
    write.csv(param_df, "param_df.csv", row.names = F)

    # Organise treatment timings for plotting:

    treat_ons <- rep(treat_ons, each = 2)
    treat_offs <- rep(treat_offs, each = 2)
    treat_offs <- c(0, treat_offs)
    
    length(treat_ons) <- min(length(treat_ons), length(treat_offs))
    length(treat_offs) <- min(length(treat_ons), length(treat_offs))
    
    treat_df <- data.frame("start" = treat_offs, "end" = treat_ons)
    treat_df["treat"] <- rep(c("off", "on"), length.out = nrow(treat_df))
    
    # Organise dataframes for plotting:

    lin_sol_dfl <- lin_sol_df %>%
    pivot_longer(cols = matches("n\\D{1}"),
                 names_to = "pheno",
                 values_to = "n")
  
    lin_sol_dfl$pheno <- factor(lin_sol_dfl$pheno, levels = c("nR", "nS", "nE"))
    
    lin_pop_df["rep"] <- rep(unique(lin_sol_df$rep), each = nrow(lin_pop_df)/4)
    
    # Remove POT samples for plotting subsetted lineages and pivot longer:

    sub_lin_df <- subset(sub_lin_df, rep != "POT")

    sub_lin_dfl <- pivot_longer(sub_lin_df,
                            cols = matches("n\\D{1}$"),
                            names_to = "pheno",
                            values_to = "n") %>% 
    mutate(pheno = str_remove(pheno, ".*_")) %>%
    dplyr::select(c("rep", "bc", "t", "N")) %>%
    unique()
    
    sub_lin_dfl$bc <- as.factor(sub_lin_dfl$bc)

    sub_pop_df <- sub_lin_dfl[, c("bc", "N", "t", "rep")] %>% unique()
    colnames(sub_pop_df) <- c("Identity", "Population", "Time", "rep")
    
    # Expand to fill missing timepoints:
    
    sub_pop_df_exp <- plyr::ddply(sub_pop_df, .(rep), function(pop_df){
      if(dim(pop_df)[1] != length(unique(pop_df$Identity)) * length(unique(pop_df$Time))) {
        added_rows <- expand.grid(Identity = unique(pop_df$Identity), Time = unique(pop_df$Time))
        added_props <- group_by(pop_df, Identity) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          dplyr::select(-one_of("Time", "Population"))
        added_rows <- merge(added_rows, added_props, all = TRUE)
        pop_df <- merge(added_rows, pop_df, all = TRUE)
        pop_df[is.na(pop_df$Population), "Population"] <- 0
        pop_df <- arrange_(pop_df, ~Time)
        warning("missing population sizes replaced by zeroes")
        return(pop_df)
      }
    })
    
    sub_pop_rf_df_exp <- sub_pop_df_exp %>% 
      dplyr::group_by(Time, rep) %>%
      dplyr::mutate(rf = Population/sum(Population))

    
    # Create a custom barcode lineage palette, ordered by most common: 

    bc_pal_df <- data.frame(bc = unique((sub_lin_df %>% arrange(desc(N)))$bc))
    cust_pal <- iwanthue(n = nrow(bc_pal_df), random = T)
    rep_pal <- colorRampPalette(brewer.pal("Paired", n = 12))(12)
    cust_pal[1:12] <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
                        "#6A3D9A", "#A6CEE3", "#FFFF99", "#B15928",
                        "#B2DF8A", "#FDBF6F", "#FB9A99", "#CAB2D6")
  
    names(cust_pal) <- bc_pal_df$bc

    # Calculate observed qD statistics

    qD_gt_P1_df <- qD_df_fxn(lin_lin_df, P="P1", q=2)
    qD_gt_P1_df["P"] <- "1"
    qD_gt_P2_df <- qD_df_fxn(lin_lin_df, P="P2", q=2)
    qD_gt_P2_df["P"] <- "2"
    
    qD_gt_df <- bind_rows(qD_gt_P1_df, qD_gt_P2_df)
    
    qD_gt_df["rep"] <- sub("(DT\\d{1})_P\\d{1}", replacement = "\\1",
                            qD_gt_df$samp)
    
    colnames(qD_gt_df) <- sub("qD", "gt_qD", colnames(qD_gt_df))

    # Lineage population plot:
    ##########################

    lin_traj_plot <- ggplot() + 
    geom_rect(data = treat_df, 
              aes(ymin = 1.0, ymax = 1e+08, 
                  xmin = start, xmax = end,
                  fill = treat), alpha = 0.2, colour = "black") +
    scale_fill_manual(values = c("palegreen1", "red"), name = "Treatment") + 
    geom_point(data = sub_lin_dfl, aes(x = t, y = N, colour = bc,
                                       group = bc),
               size = 0.4) +
    geom_line(data = unique(lin_sol_dfl %>% 
                              dplyr::select(t, N, rep) %>%
                              dplyr::mutate(rep = paste0("DT", rep))),
              aes(x = t, y = N),
              size = 1.2, linetype = "dotted") +
    geom_point(data = lin_pop_df %>%
                 subset(t != 4.0) %>%
                 dplyr::mutate(rep = paste0("DT", rep)),
               aes(x = t, y = u),
               size = 3.0) +
    scale_colour_manual(values = cust_pal) +
    xlab("Time (Days)") + 
    ylab("log(Population)\n") +
    scale_x_continuous(limits = c(0, 80.0)) +
    scale_y_log10(limits = c(1, 1e+08)) + 
    facet_wrap(~rep, ncol = 1) +
    theme_BARCODE(text_size=18) + 
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank()
    )

    # Lineage frequency plot: 
    #########################

    lin_freq_plot <- ggplot(data = sub_pop_rf_df_exp, aes(x = Time, y = rf,
                                                          fill = as.factor(Identity))) +
        geom_area() +
        scale_x_continuous(limits = c(0, 80.0),
                            name = "Time (Days)") +
        ylab("Lineage\nRelative Frequency\n") +
        scale_fill_manual(values = cust_pal) +
        facet_wrap(~rep, ncol = 1) +
        theme_BARCODE(text_size=18) + 
        theme(legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_blank()
        )

    # qD statistic plot
    ###################
    
    qD_plot <- ggplot(data = qD_gt_df, aes(x = gt_qD, y = gt_qD_diss, fill = gt_qD_diss)) + 
        geom_point(shape = 21, colour = "black", alpha = 0.8, size = 5) + 
        scale_x_log10(limits = c(1e+00, 1e+06),
                    breaks = c(1e+01, 1e+03, 1e+05),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) + 
        scale_y_continuous(limits = c(0.5, 4.5),
                            breaks = 1:4,
                            labels = 1:4) + 
        xlab("Lineage\nDiversity\n") + 
        ylab("Lineage Diversity Dissimilarity") +
        scale_fill_viridis_c(limits = c(0.0, 4.0),
                            name=  "Diversity\nDissimilarity") +
        facet_grid(rep~P) + 
        theme_BARCODE(text_size=18) + 
        theme(panel.grid.major = element_line(size = 0.5, colour = "grey90"),
            panel.border = element_rect(colour = "black", fill=NA, size=1.0),
            axis.text.x = element_text(size = 12))

    # Combined Plot
    ###############

    combined_plot <- cowplot::plot_grid(lin_traj_plot,
                                        lin_freq_plot,
                                        qD_plot,
                                        nrow = 1, rel_widths = c(1.0,1.0,0.8),
                                        axis = "bt",
                                        align = "hv")

    # Save the plot
    ###############

    # Change ggsave function so bg = white
    ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

    ggsave("comb_sub_lin_traj_plot.pdf",
           plot = combined_plot,
           width = 19.0, height = 7.0)

    ggsave("comb_sub_lin_traj_plot.jpg",
           plot = combined_plot,
           width = 19.0, height = 7.0)

"""


# Print final time
println(now())

################################################################################

