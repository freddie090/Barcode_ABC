
using CSV
using RCall
using DataFrames

# Simulates the agent-based lineage model given the most recent pyABC 
# generation generated using the change in cell population sizes. 
# The Simulation type (Model A - Unidirectional Transitions, 
#                      Model B - Bidirectional Transitions, 
#                      Model C - Escape Transitions), 
# simulation output directory (e.g. 'sim_1') and job_id (one for every 
# draw from the first step of the fitting process) are loaded from respective 
# command line arguments.

length(ARGS) == 5 || error("Should be 5 command line arguments.")

abc_mod_dir = ARGS[1]
abc_out_dir = ARGS[2]
abc_model = ARGS[3]
out_dir_name = ARGS[4]
job_id = Int64(eval(Meta.parse(ARGS[5])))
sim_type = ARGS[6]
param_df_name = ARGS[7]

abc_model ∈ ["A", "B", "C"] || error("The 'abc_model' variable is not one of the three permitted values: A, B or C.")

# Load in functions for the lineage ('Barcode') simulations. 
include(string(abc_mod_dir, "/ABC_Barcode_KMC_Model.jl"))

# Print current time to output
println(now())
# Check number of threads being used.
println(Base.Threads.nthreads())

cd(abc_out_dir)
cd(out_dir_name)

# Load in the final pyABC generation:
param_df = DataFrame(CSV.File(param_df_name))

# Use the job id to choose which row to use
param_df_samp = param_df[job_id,:]

# Extract parameters
rho = 10^(-param_df_samp[:rho])
mu = 10^(-param_df_samp[:mu])
Dc = param_df_samp[:Dc]
if "delt" in names(param_df_samp)
    del = 10^(-param_df_samp[:delt])
else
    del = param_df_samp[:del]
end
k = param_df_samp[:k]
psi = param_df_samp[:psi]

if abc_model == "A"
    sig = 0.0
    al = 0.0
end

if abc_model == "B"
    sig = 10^(-param_df_samp[:sig])
    al = 0.0
end

if abc_model == "C"
    sig = 10^(-param_df_samp[:sig])
    al = 10^(-param_df_samp[:al])
end


################################################################################

########################################################
# Set the fixed parameters for the simualted experiment:
######################################################## 

# Number of uniquely barcoded cells when the experiment :
n0 = Int64(1e+04);

# Estimates of the average, sensitive population's birth and death rates:
b = 0.893; d = 0.200;

# The measurement error of the population size readings (where
# observed N ~ Normal(mean = true N, sd = true N * phi)):
phi = 0.10;

# Time when the experiment begins and maximum time of the experiment (in days):
t0 = 0.0; tmax = 60.0;

# Maximum number of cells before a replicate is harvested, and an estimate of
# the carrying capacity:begins
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

################################################################################

lin_sim_out = exp_pulse_treat_kmc_full_sim_fxn(n0, b, d, rho, mu, sig, del, al,
                                                Dc, k, psi, 
                                                t0, t_exp, tmax, t_Pass, 
                                                Nmax, Cc,
                                                treat_ons, treat_offs, 
                                                dt_save_at, t_keep,
                                                n_Pass=n_Pass, R_real="l")

# Barcode lineage distributions 
lin_lin_df = lin_sim_out["lin_df"]
lin_pop_df = DataFrame(t = lin_sim_out["t"], u = lin_sim_out["u"])
lin_sol_df = lin_sim_out["sol_df"]

@rput lin_lin_df;
@rput lin_pop_df;
@rput lin_sol_df;
@rput job_id; @rput sim_type;
                        
R"""

    # Save dataframes
    write.csv(lin_lin_df, paste0("KMC_", sim_type, "_lin_df_", job_id, ".csv"), 
            row.names = F)
    write.csv(lin_pop_df, paste0("KMC_", sim_type, "_pop_df_", job_id, ".csv"), 
            row.names = F)
    write.csv(lin_sol_df, paste0("KMC_", sim_type, "_sol_df_", job_id, ".csv"), 
            row.names = F)

"""

# Print final time
println(now())

################################################################################

end