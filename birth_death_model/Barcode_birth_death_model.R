
# Copyright 2025 Cancer Research Technology and The Institute of Cancer Research.
#
# Licensed under a software academic use license provided with this software package (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at: https://github.com/freddie090/Barcode_ABC
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

################################################################################

# Barcode Sampling Noise Birth Death Model

library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set the working directory (update as appropriate).
setwd("/path/to/working/directory")

# Load custom functions for model fitting.
source("/path/to/Barcode_birth_death_custom_functions.R")

# Load count frequency tables for analysis.
setwd("Count_Tables")

bc_dfs <- lapply(grep(".*_count_freq_table.csv", list.files(), value = TRUE), read.csv, stringsAsFactors = FALSE)

# Format the data.
bc_dfs <- lapply(bc_dfs, function(x) {
  colnames(x) <- c("count", "freq")
  return(x)
})

Js <- lapply(bc_dfs, function(bc_df) {
  sum(bc_df$count * bc_df$freq)
})

# Define parameters.
t <- 7.0       # Duration of experiment (in days).
K <- 10^6      # Number sampled during binomial sampling step.
I <- length(bc_dfs)  # Number of replicate samples.
nmax <- 1000   # Maximum p(k | n) to calculate.
Nt <- 80 * (10^6) # Number of cells post-expansion.
N0 <- 10^6    # Estimated number of uniquely barcoded cells at t=0.

# Rescale and compress counts for each replicate.
bc_dfs <- lapply(bc_dfs, function(bc_df) {
  J <- sum(bc_df$count * bc_df$freq)
  KJ_re_scale <- K / J
  
  bc_df$count <- round(bc_df$count * KJ_re_scale)
  bc_df <- subset(bc_df, count > 0)
  bc_df <- aggregate(freq ~ count, data = bc_df, FUN = sum)

  bc_df <- data.frame(
    count = c(0, bc_df$count), freq = c((N0 - sum(bc_df$freq)), bc_df$freq)
  )
  return(bc_df)
})

# Prepare data for Stan model.
Ncv <- unlist(lapply(bc_dfs, nrow))
csv <- unlist(lapply(bc_dfs, function(x) { x$count }))
fsv <- unlist(lapply(bc_dfs, function(x) { x$freq }))
Nc <- length(csv)
Nvec <- 1:Nc

# Compute variance matrix.
kmax <- 1000
jvar_mat <- matrix(nrow = nmax + 1, ncol = I)

for (i in seq_along(Js)) {
  for (x in 1:(nmax + 1)) {
    jvar_mat[x, i] <- get_normalised_variance_with_zeros((x - 1), Nt, K, Js[[i]], kmax)
  }
}

# Prepare data list for Stan.
bc_df_dat <- list(N0 = N0, Nc = Nc, Nvec = Nvec, t = t, K = K, I = I, 
                  nmax = nmax, Nt = Nt, Ncv = Ncv, csv = csv, fsv = fsv, 
                  jvar_mat = jvar_mat)

# Fit the model.
bc_df_fit <- stan(file = "Barcode_birth_death_model.stan", 
                  data = bc_df_dat, chains = 4, iter = 4000)

# Save the output.
saveRDS(bc_df_fit, file = "Barcode_birth_death_model_fit.rds")
