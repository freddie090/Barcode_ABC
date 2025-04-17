
using RCall
using DataFrames

# Simulates synthetic lineage sizes and population size changes through 
# treatment given known parameter values. Uses both models to simulate
# total population size changes through treatment (hybrid ODE/stochastic jump 
# and agent-based kinetic monte carlo (KMC)) or lineage size Distributions
# (just agent-based kinetic monte carlo (KMC)).
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

################################################################################

###############
# Run the Model
###############

# To make sure that the (cheaper) phenotypic compartment model is accurately
# recovering the (more memory intensive) agent-based lineage model, we run 
# both in parallel and save the outputs from each alongside a plot of the 
# simulated population trajectories.
# (n.b. the simulations are stochastic so the outputs might not look
# identical)

# Total population ODE version

ode_sim_out = exp_pulse_treat_full_sol_fxn(n0, b, d, rho, mu, sig, del, al, 
                                           Dc, k, psi,
                                           t0, t_exp, tmax, t_Pass, 
                                           Nmax, Cc,
                                           treat_ons, treat_offs, 
                                           t_keep, Nswitch, 
                                           save_at=0.2, n_Pass=n_Pass)

# Add measurement noise to the population sizes
ode_sim_out["u"] = model_meas_noise(ode_sim_out["u"], phi)

# Now run the analagous model, but the barcode lineage version (simulates
# each barcoded lineage independently, opposed to just sensitive and
# resistant populations) - also saves full solution for N. 

lin_sim_out = exp_pulse_treat_kmc_full_sim_fxn(n0, b, d, rho, mu, sig, del, al, 
                                                Dc, k, psi, 
                                                t0, t_exp, tmax, t_Pass, 
                                                Nmax, Cc,
                                                treat_ons, treat_offs, 
                                                dt_save_at, t_keep,
                                                n_Pass=n_Pass, R_real="l")


lin_sim_out["u"] = model_meas_noise(lin_sim_out["u"], phi)


# Now, save the outputs for:
#  i) Population lineage fitting via ABC
# ii) Lineage distribution comparisons via qD statistics

# ODE observed populations and parameters
ode_pop_df = DataFrame(t = ode_sim_out["t"], u = ode_sim_out["u"])

# ODE full solution population trajectories
ode_sol_df = ode_sim_out["sol_df"]

# Barcode lineage observed populations and parameters
lin_pop_df = DataFrame(t = lin_sim_out["t"], u = lin_sim_out["u"])

# Barcode lineage full solution population trajectories
lin_sol_df = lin_sim_out["sol_df"]

# Barcode lineage distribution 
lin_lin_df = lin_sim_out["lin_df"]

# This simulation's parameters 
param_df = DataFrame(abc_model = abc_model, 
                     rho = rho, mu = mu, sig = sig, al = al, del = del,                        
                     Dc = Dc, k = k, psi = psi, phi = phi,
                     n0 = n0, b = b, d = d, Nmax = Nmax, Cc = Cc)

# Use RCall to save dataframes and generate plots.
@rput ode_pop_df; @rput ode_sol_df;
@rput lin_pop_df; @rput lin_lin_df;
@rput lin_sol_df;
@rput param_df;
@rput treat_ons; @rput treat_offs; @rput tmax;
@rput abc_model; @rput sim_id;

# Save the simulation output dataframes and phenotype population trajectory
# plots for both model types.

R"""

    library(ggplot2)

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
    write.csv(ode_pop_df, "ode_pop_df.csv", row.names = F)
    write.csv(ode_sol_df, "ode_sol_df.csv", row.names = F)
    write.csv(lin_pop_df, "lin_pop_df.csv", row.names = F)
    write.csv(lin_lin_df, "lin_lin_df.csv", row.names = F)
    write.csv(lin_sol_df, "lin_sol_df.csv", row.names = F)
    write.csv(param_df, "param_df.csv", row.names = F)

    # Save the population trajectories for both the ode and kmc simulations: 

    treat_ons <- rep(treat_ons, each = 2)
    treat_offs <- rep(treat_offs, each = 2)
    treat_offs <- c(0, treat_offs)

    length(treat_ons) <- min(length(treat_ons), length(treat_offs))
    length(treat_offs) <- min(length(treat_ons), length(treat_offs))

    treat_df <- data.frame("start" = treat_offs, "end" = treat_ons)
    treat_df["treat"] <- rep(c("off", "on"), length.out = nrow(treat_df))

    theme_BARCODE <- function(){
    theme_minimal() +
        theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    }

    lin_plot <- ggplot() + 
        geom_rect(data = treat_df, 
                aes(ymin = 1.0, ymax = 1e+09, 
                    xmin = start, xmax = end,
                    fill = treat), alpha = 0.2, colour = "black") +
        scale_fill_manual(values = c("palegreen1", "red"), name = "Treatment") + 
        geom_line(data = lin_sol_df, aes(x = t, y = N, group = rep),
                alpha = 0.2, size = 0.8) +
        geom_point(data = lin_pop_df, 
                aes(x = t, y = u), size = 3) + 
        scale_y_log10(limits = c(1, 1e+09)) + 
        annotation_logticks(base = 10, sides = "l") +
        scale_x_continuous(limits = c(0, tmax)) + 
        xlab("Time (days)") +
        ylab("N\n") +
        theme_BARCODE() + 
        ggtitle("KMC Simulation")

        ode_plot <- ggplot() + 
        geom_rect(data = treat_df, 
                aes(ymin = 1.0, ymax = 1e+09, 
                    xmin = start, xmax = end,
                    fill = treat), alpha = 0.2, colour = "black") +
        scale_fill_manual(values = c("palegreen1", "red"), name = "Treatment") + 
        geom_line(data = ode_sol_df, aes(x = t, y = N, group = rep),
                alpha = 0.2, size = 0.8) +
        geom_point(data = ode_pop_df, 
                aes(x = t, y = u), size = 3) + 
        scale_y_log10(limits = c(1, 1e+09)) + 
        annotation_logticks(base = 10, sides = "l") +
        scale_x_continuous(limits = c(0, tmax)) + 
        xlab("Time (days)") +
        ylab("N\n") +
        theme_BARCODE() + 
        ggtitle("ODE Simulation")

    # Change ggsave function so bg = white
    ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

    ggsave("lin_pop_traj_plot.jpg", lin_plot, height = 4, width = 12)
    ggsave("ode_pop_traj_plot.jpg", ode_plot, height = 4, width = 12)

"""

# Print final time
println(now())

################################################################################

