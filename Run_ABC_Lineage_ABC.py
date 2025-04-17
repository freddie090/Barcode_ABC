
# Run the first step of the ABC parameter inference step that uses the 
# total cell population size changes through treatment.

import os 
import tempfile
import sys
from datetime import timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random

# log10 fxn for plotting
import math
log10v = np.vectorize(math.log10)

import pyabc
from pyabc import ABCSMC, RV, Distribution, SingleCoreSampler
from pyabc.external.julia import Julia

pyabc.settings.set_figure_params('pyabc')  # for beautified plots

def str_to_bool(s: str) -> bool:
    return s.strip().lower() in ['true', '1', 'yes', 'y']

# Load in command line arguments:
abc_mod_dir = sys.argv[1]
abc_out_dir = sys.argv[2]
abc_model = sys.argv[3]
abc_fit_dir = sys.argv[4]
is_simulated_data = str_to_bool(sys.argv[5])

if abc_model not in ["A", "B", "C"]: 
    raise ValueError("The 'abc_model' variable is not one of the three permitted values: A, B or C.")

os.chdir(abc_out_dir)
os.chdir(abc_fit_dir)

# Read in the output dataframes
lin_pop_df = pd.read_csv('lin_pop_df.csv')

# Custom log function for reading in parameters to account for models where
# parameters have been set to 0.0: 

def custom_log(x):
    if x > 0:
        return -math.log10(x)
    elif x == 0:
        return -1
    else:
        raise ValueError("Input must be non-negative")

################################################################################

###########
# Setup ABC
###########

# if is_simulated_data:
param_df = pd.read_csv('param_df.csv')

rho = param_df.iloc[0]['rho']
mu = param_df.iloc[0]['mu']
delt = param_df.iloc[0]['del']
Dc = param_df.iloc[0]['Dc']
k = param_df.iloc[0]['k']
psi = param_df.iloc[0]['psi']

if abc_model == "C":
    al = param_df.iloc[0]['al']
    sig = param_df.iloc[0]['sig']
    gt_par = {"rho" : custom_log(rho), "mu": custom_log(mu), 
                "sig": custom_log(sig), "al" : custom_log(al),
                "Dc": Dc, "delt" : delt, "k" : k, "psi" : psi}
elif abc_model == "B":
    al = 0.0
    sig = param_df.iloc[0]['sig']
    gt_par = {"rho" : custom_log(rho), "mu": custom_log(mu), 
                "sig": custom_log(sig),
                "Dc": Dc, "delt" : delt, "k" : k, "psi" : psi}
else:
    al = 0.0
    sig = 0.0
    gt_par = {"rho" : custom_log(rho), "mu": custom_log(mu), 
        "Dc": Dc, "delt" : delt, "k" : k, "psi" : psi}


# Load the julia phenotypic compartment model used for ABC: 
jl = Julia(module_name="ABC_Population_Model", 
        source_file=os.path.join(abc_mod_dir, "ABC_Population_Model.jl"))

# Load module model
model = jl.model()

# Load module distance function
distance = jl.distance()


# precompile model function
model(gt_par)
# and distance function
distance(model(gt_par), model(gt_par))

# Load in the observation from the lineage population trajectory dataframe
obs = {"t" : lin_pop_df['t'].values, "u" : lin_pop_df['u'].values}
distance(obs, model(gt_par))

# Set prior distributions:

prior = Distribution(
    rho = RV("uniform", 0, 7),
    mu = RV("uniform", 0, 9),
    Dc = RV("gamma", 3, scale=1/1),
    delt = RV("uniform", 0, 3),
    k = RV("gamma", 2, scale=1/1),
    psi = RV("beta", 3, 1)
)

if abc_model == "B":
    prior["sig"] = RV("uniform", 0, 9)

if abc_model == "C":
    prior["sig"] = RV("uniform", 0, 9)
    prior["al"] = RV("uniform", 0, 9)

################################################################################

################
# Run ABC Step 1
################

# Number of populations to retain per generation
pop_size = 100

# Create pyabc object
abc = ABCSMC(
    model,
    prior,
    distance,
    sampler=SingleCoreSampler(),
    population_size = pop_size
)

# Set the .db path
db_path = os.path.join(abc_out_dir, abc_fit_dir, 
                    f"ABC_Pop_Inf_out.db")

# Run the inference
history = abc.new("sqlite:///" + db_path, obs)
history = abc.run(max_nr_populations=6,
                  minimum_epsilon=1e-02,
                  max_walltime=timedelta(hours = 94))

################################################################################

#################
# Save ABC Output
#################

# Save the final generation as a dataframe

h = pyabc.storage.History("sqlite:///"+f"ABC_Pop_Inf_out.db")

df = h.get_distribution()[0]

filename = f"ABC_Pop_final_gen_dist.csv"

df.to_csv(filename)

################################################################################

