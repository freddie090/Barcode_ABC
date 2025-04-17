
################################################################################

Code to run the simulations and two-step ABC inference from the manuscript: 

"Quantitative measurement of phenotype dynamics during cancer drug resistance 
evolution using genetic barcoding"

Correspondence: freddie.whiting@icr.ac.uk

################################################################################

#####################
# Script Descriptions
#####################

The directory contains the following scripts:

1. ABC_Population_Model.jl [julia]

Combines a deterministic ODE with stochastic jump process to model the 
evolution of 2/3 (depending on model A,B,C) phenotypic compartments through
treatment. This module is loaded as part of the first ABC step (via pyABC)
in the python script 'Run_ABC_Lineage_ABC.py'. Details of how the parameters
that control the evolution of the phenotypes are described in full in the 
manuscript.


2. ABC_Population_Model_Fxns.jl [julia]

Contains the functions that make up 1., but in a format so they can be used 
for the generation of synthetic data.


3. ABC_Barcode_KMC_Model.jl [julia]

The kinetic monte-carlo (KMC) agent-based model that simulates the analogous  
models to those described in 1., but now keeps track of each individual cell 
lineage. Returns the cell lineage distributions for each 'experimental passage',
per replicate.


4. Run_ABC_Lineage_Sim.jl [julia]

Simulates synthetic data that mirrors data captured during a lineage tracing
experimental resistance evolution setup. Namely, total population size changes
through treatment at specified times (denoted by t_keep) including passage 
times (t_Pass) and the time when the experiment ends (tmax). Premature 
passaging/experiment ends occur if the populations exceed the maximum population
size (Nmax) or experience extinction. Parameters that control the evolution of 
the cell phenotypes are loaded in from a parameter table csv.


5. Rub_ABC_Lineage_Lin_Track.jl [julia]

Simulates the same data as in 'Run_ABC_Lineage_Sim.jl', but now also keeps
track of the most frequent x lineages regularly throughout the simulation. 
This enables plotting of the full lineage trajectories through treatment 
(as in Fig.2 of the manuscript). These plots are automatically generated and 
saved after the simulation has completed.


6. Run_ABC_Lineage_ABC.py [python]

Performs the first stage of the ABC inference - performs ABC-SMC using just the 
total cell population size observations through treatment. 
Saves the output of the fitting process as a .csv of the final generation.


7. Run_ABC_Lineage_Fit.jl [julia]

Takes the parameter distribution from the first step of the ABC fitting 
process (only uses the total population sizes) and performs additional 
simulations using the agent-based model (KMC) version. 'job_id' (loaded
from the command line) dictates which line of the posterior is used for any
given simulation's parameters. Can also be run with the final posterior
distribution to generate posterior predictive simulations (PPS). 


8. Run_ABC_Posterior_Plot_Inference.R [R]

Takes the simulated outputs from 7. and produces summary plots and the final
posterior distribution by retaining the parameter values that generated
lineage diversity statistics closest to the observed values. Can also be set
to only generate plots if using with the posterior predictive simulations 
outputs.


9. Barcode_ABC_Data_Functions.R [R]
10. Barcode_ABC_Plotting_Functions.R [R]

Functions for manipulating and plotting the output data in R.


################################################################################

#####################################
Installation and System Requirements:
#####################################

The computational code for this project requires the following software 
environments and dependencies across Python, Julia, and R. The code assumes 
compatibility with the following versions and packages:

Python Requirements:

Python Version: 3.1 or later.
Python Libraries:
matplotlib: For plotting and visualisation (tested with version 3.7+).
numpy: For numerical computations (version 1.24+).
pandas: For data manipulation (version 2.0+).
pyabc: For Approximate Bayesian Computation (version 0.16+).
math: Standard library for mathematical operations.
os, tempfile, sys, random: Standard library for file handling and system 
operations.
pyabc.external.julia: Interface to Julia (compatible with PyJulia).
Note: The code requires the PyJulia interface, and Julia must be installed and 
accessible from the Python environment.


Julia Requirements:

Julia Version: 1.7 or later.
Julia Packages:
DifferentialEquations.jl: For solving ODEs and SDEs (version 7.6+).
JumpProcesses.jl: For stochastic process simulation.
Distributions.jl: Statistical distributions (version 0.25+).
DataFrames.jl: Data manipulation (version 1.3+).
CSV.jl: CSV file input/output (version 0.10+).
StatsBase.jl: Basic statistical functions.
FreqTables.jl: Frequency table computation.
RCall.jl: Interface for calling R from Julia.
Dates: Standard library for date and time handling.
Suppressor.jl: For suppressing output in scripts.


R Requirements:

R Version: 4.2 or later.
R Packages:
ggplot2: For data visualisation (version 3.4+).
dplyr, tidyr, plyr: Data manipulation packages.
cowplot: For combining plots.
RColorBrewer: Colour palettes for visualisation.
scales: Scaling functions for visualisation.
reshape2: For data reshaping.
stringr: For string manipulation.
purrr: Functional programming tools.
Cairo: For improved graphics device support.


Dependencies Installation:

Python packages can be installed via pip or conda.

Julia packages can be added via Julia's package manager 
(e.g., Pkg.add("DifferentialEquations")).

R packages can be installed using install.packages().

Cross-language Integration:
The software requires cross-language integration via PyJulia (Python to Julia) 
and RCall (Julia to R). This necessitates proper installation and configuration 
of all three languages in the same environment, with R accessible from Julia and
Julia accessible from Python.

################################################################################

########
# Usage:
########

To run the simulation and inference framework using the cell numbers described 
in the paper is computationally expensive and should therefore be run in 
parallel on a HPC.

To use the framework on a new dataset, fixed parameters that describe the 
experimental design used for data generation must be modified in the following
scripts:
- `ABC_Population_Model.jl`
- `Run_ABC_Lineage_Sim.jl`
- `Run_ABC_Lineage_Fit.jl`

For testing purposes, synthetic data with known parameters, reduced cell 
numbers, and fewer simulation iterations can be run locally. An example using 
'Model A  Unidirectional Transitions' from the paper is provided below.

A parameter table for generating example data for 'Model A' is included 
alongside the scripts:
- 'ABC_Model_A_Param_Table.csv'

```shell

# Define paths and parameters
abc_script_dir = "...path/to/abc_script_dir"
abc_mod_dir = "...path/to/model_dir"
abc_out_dir = "...path/to/abc_out_dir"
abc_model = "A"
param_table_path = "...path/to/ABC_Model_A_Param_Table.csv"
sim_id = "1"
abc_fit_dir = "Model_A_Outputs/sim_1"
is_simulated_data = "True"
post_n = "20"


# Step 1: Run the simulation to generate synthetic data

julia "${abc_script_dir}/Run_ABC_Lineage_Sim.jl" \
  $abc_mod_dir $abc_out_dir $abc_model $param_table_path $sim_id


# Step 2: Run the first step of the inference using the population size changes

python-jl "${abc_script_dir}/Run_ABC_Lineage_ABC.py" \
  $abc_mod_dir $abc_out_dir $abc_model $abc_fit_dir $is_simulated_data


# Step 3: Simulate the agent-based lineage simulation using parameters from 2.

sim_type = "Sim"
param_df_name = "ABC_Pop_final_gen_dist.csv"

for job_id in {1..100}; do 
    julia "${abc_script_dir}/Run_Lineage_Fit.jl" \
      $abc_mod_dir $abc_out_dir $abc_model $abc_fit_dir $job_id $sim_type \
      $param_df_name
done


# Step 4: Calculate the posterior distribution and generate plots

calculate_posterior = "T"

Rscript "${abc_script_dir}/Run_ABC_Posterior_Plot_Inference.R" \
  $abc_script_dir $abc_mod_dir $abc_out_dir $abc_model $abc_fit_dir \
  $output_type $is_simulated_data $param_df_name $calculate_posterior $post_n


# Step 5: Run the posterior predictive simulations (PPS)

sim_type = "PPS"
param_df_name = "ABC_post_dist.csv"

for job_id in {1..20}; do 
    julia "${abc_script_dir}/Run_Lineage_Fit.jl" \
      $abc_mod_dir $abc_out_dir $abc_model $abc_fit_dir $job_id $sim_type \ 
      $param_df_name
done


# Step 6: Generate plots of the posterior predictive simulations (PPS)

calculate_posterior = "F"

Rscript "${abc_script_dir}/Run_ABC_Posterior_Plot_Inference.R" \
  $abc_script_dir $abc_mod_dir $abc_out_dir $abc_model $abc_fit_dir \
  $output_type $is_simulated_data $param_df_name $calculate_posterior

```

The example will produce a posterior distribution for the model's parameters and
plots to assess the quality of the fit by showing the inferred vs observed 
population size changes and lineage diversity statistics.

################################################################################