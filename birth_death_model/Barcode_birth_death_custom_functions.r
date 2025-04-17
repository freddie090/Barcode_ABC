
# Copyright 2025 Cancer Research Technology and The Institute of Cancer Research.
#
# Licensed under a software academic use license provided with this software package (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at: https://github.com/freddie090/Barcode_ABC
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

################################################################################

#' Expectation of Zero-Truncated Poisson
#'
#' Calculates the conditional expectation of a zero-truncated Poisson distribution.
#'
#' @param k Observed count.
#' @param K Total count.
#' @param J Normalisation factor.
#' @return Conditional expectation of the zero-truncated Poisson.
#' @examples
#' expectation_zt_poisson(5, 10, 20)
expectation_zt_poisson <- function(k, K, J){
  mu <- (k / K) * J
  cE <- mu / (1 - exp(-mu))
  return(cE)
}

#' Variance of Zero-Truncated Poisson
#'
#' Calculates the conditional variance of a zero-truncated Poisson distribution.
#'
#' @param k Observed count.
#' @param K Total count.
#' @param J Normalisation factor.
#' @return Conditional variance of the zero-truncated Poisson.
#' @examples
#' variance_zt_poisson(5, 10, 20)
variance_zt_poisson <- function(k, K, J){
  mu <- (k / K) * J
  mu_ztP <- expectation_zt_poisson(k, K, J)
  p_j0 <- dpois(0, mu)
  cV <- mu / (1 - p_j0) - (p_j0 * mu_ztP^2)
  return(cV)
}

#' Expectation of Zero-Truncated Negative Binomial
#'
#' Calculates the conditional expectation of a zero-truncated negative binomial distribution.
#'
#' @param k Observed count.
#' @param K Total count.
#' @param J Normalisation factor.
#' @param phi Dispersion parameter.
#' @return Conditional expectation of the zero-truncated negative binomial.
#' @examples
#' expectation_zt_negbinom(5, 10, 20, 0.5)
expectation_zt_negbinom <- function(k, K, J, phi){
  mu <- (k / K) * J
  cE <- mu / (1 - (1 + (mu / phi))^(-phi))
  return(cE)
}

#' Variance of Zero-Truncated Negative Binomial
#'
#' Calculates the conditional variance of a zero-truncated negative binomial distribution.
#'
#' @param k Observed count.
#' @param K Total count.
#' @param J Normalisation factor.
#' @param phi Dispersion parameter.
#' @return Conditional variance of the zero-truncated negative binomial.
#' @examples
#' variance_zt_negbinom(5, 10, 20, 0.5)
variance_zt_negbinom <- function(k, K, J, phi){
  mu <- (k / K) * J
  mu_ztODP <- expectation_zt_negbinom(k, K, J, phi)
  p_j0 <- dnbinom(0, mu = mu, size = phi)
  cV <- (mu * (1 + (mu / phi))) / (1 - p_j0) - (p_j0 * mu_ztODP^2)
  return(cV)
}

#' Compound Variance of Zero-Truncated Poisson
#'
#' Calculates the compound variance by marginalising over all possible values of k.
#'
#' @param n Observed lineage size.
#' @param Nt Total population size.
#' @param K Total count.
#' @param J Normalisation factor.
#' @param kmax Maximum value of k.
#' @return Compound variance of the zero-truncated Poisson.
#' @examples
#' compound_variance_zt_poisson(50, 1000, 10, 20, 100)
compound_variance_zt_poisson <- function(n, Nt, K, J, kmax){
  p_ks <- dpois(0:kmax, (K * (n / Nt)))
  vec1 <- sapply(0:kmax, variance_zt_poisson, K = K, J = J) * p_ks
  vec1[is.na(vec1)] <- 0
  vec2 <- sapply(0:kmax, expectation_zt_poisson, K = K, J = J)^2 * p_ks
  vec2[is.na(vec2)] <- 0
  vec3 <- sapply(0:kmax, expectation_zt_poisson, K = K, J = J) * p_ks
  vec3[is.na(vec3)] <- 0
  comp_var <- sum(vec1) + (sum(vec2) - sum(vec3)^2)
  return(comp_var)
}

#' Compound Variance of Zero-Truncated Negative Binomial
#'
#' Calculates the compound variance by marginalising over all possible values of k.
#'
#' @param n Observed lineage size.
#' @param Nt Total population size.
#' @param K Total count.
#' @param J Normalisation factor.
#' @param kmax Maximum value of k.
#' @param phi Dispersion parameter.
#' @return Compound variance of the zero-truncated negative binomial.
#' @examples
#' compound_variance_zt_negbinom(50, 1000, 10, 20, 100, 0.5)
compound_variance_zt_negbinom <- function(n, Nt, K, J, kmax, phi){
  p_ks <- dpois(0:kmax, (K * (n / Nt)))
  vec1 <- sapply(0:kmax, variance_zt_negbinom, K = K, J = J, phi = phi) * p_ks
  vec1[is.na(vec1)] <- 0
  vec2 <- sapply(0:kmax, expectation_zt_negbinom, K = K, J = J, phi = phi)^2 * p_ks
  vec2[is.na(vec2)] <- 0
  vec3 <- sapply(0:kmax, expectation_zt_negbinom, K = K, J = J, phi = phi) * p_ks
  vec3[is.na(vec3)] <- 0
  comp_var <- sum(vec1) + (sum(vec2) - sum(vec3)^2)
  return(comp_var)
}

#' Normalised Variance for Poisson-Based Sequencing
#'
#' Pre-allocates a vector of variances for a normalised distribution assuming Poisson sequencing.
#'
#' @param n Observed lineage size.
#' @param Nt Total population size.
#' @param K Total count.
#' @param J Normalisation factor.
#' @param kmax Maximum value of k.
#' @return Vector of normalised variances.
#' @examples
#' get_normalised_variance(50, 1000, 10, 20, 100)
get_normalised_variance <- function(n, Nt, K, J, kmax){
  p_ks <- dpois(0:kmax, (K * (n / Nt)))
  vec1 <- sapply(0:kmax, variance_zt_poisson, K = K, J = J) * p_ks
  vec1[is.na(vec1)] <- 0
  vec2 <- sapply(0:kmax, expectation_zt_poisson, K = K, J = J)^2 * p_ks
  vec2[is.na(vec2)] <- 0
  vec3 <- sapply(0:kmax, expectation_zt_poisson, K = K, J = J) * p_ks
  vec3[is.na(vec3)] <- 0
  var_j <- sum(vec1) + (sum(vec2) - (sum(vec3))^2)
  var_j <- var_j * (K / J)^2
  return(var_j)
}

#' Normalised Variance Including Zero Counts
#'
#' Calculates the variance including zero counts for normalisation.
#'
#' @param n Observed lineage size.
#' @param Nt Total population size.
#' @param K Total count.
#' @param J Normalisation factor.
#' @param kmax Maximum value of k.
#' @return Vector of normalised variances including zero counts.
#' @examples
#' get_normalised_variance_with_zeros(50, 1000, 10, 20, 100)
get_normalised_variance_with_zeros <- function(n, Nt, K, J, kmax){
  p_ks <- dpois(0:kmax, (K * (n / Nt)))
  vec1 <- sapply(0:kmax, function(x) {(x / K) * J}) * p_ks
  vec1[is.na(vec1)] <- 0
  vec2 <- sapply(0:kmax, function(x) {(x / K) * J})^2 * p_ks
  vec2[is.na(vec2)] <- 0
  vec3 <- sapply(0:kmax, function(x) {(x / K) * J}) * p_ks
  vec3[is.na(vec3)] <- 0
  var_j <- sum(vec1) + (sum(vec2) - (sum(vec3))^2)
  var_j <- var_j * (K / J)^2
  return(var_j)
}

#' Probability of Lineage Size
#'
#' Computes the probability of a lineage size n given birth and death rates and time.
#'
#' @param n Lineage size.
#' @param b Birth rate.
#' @param d Death rate.
#' @param t Time.
#' @return Probability of the lineage size.
#' @examples
#' probability_lineage_size(5, 0.1, 0.05, 10)
probability_lineage_size <- function(n, b, d, t){
  alpha <- d * (exp((b - d) * t) - 1) / (b * exp((b - d) * t) - d)
  beta <- b * (exp((b - d) * t) - 1) / (b * exp((b - d) * t) - d)
  if(n == 0){
    pn <- alpha
  } else {
    pn <- (1 - alpha) * (1 - beta) * (beta^(n - 1))
  }
  return(pn)
}

#' Probability of Observing Lineage Reads
#'
#' Calculates the probability of observing lineage reads given lineage size and other parameters.
#'
#' @param n Lineage size.
#' @param x Observed reads.
#' @param b Birth rate.
#' @param d Death rate.
#' @param t Time.
#' @param K Total count.
#' @param nmax Maximum lineage size.
#' @param Nt Total population size.
#' @param variance_vector Precomputed variance vector.
#' @return Probability of observing lineage reads.
#' @examples
#' probability_lineage_reads(10, 5, 0.1, 0.05, 10, 100, 50, 1000, rep(1, 51))
probability_lineage_reads <- function(n, x, b, d, t, K, nmax, Nt, variance_vector){
  if(n == 0){
    if(x == 0){
      pxn <- 1.0
    } else {
      pxn <- 0.0
    }
  } else {
    mu <- K * (n / Nt)
    phi <- (mu^2) / (variance_vector[n + 1] - mu)
    pxn <- dnbinom(x, mu = mu, size = phi)
  }
  return(pxn)
}

#' Marginal Probability of Observing Reads
#'
#' Computes the marginal probability of observing reads by summing over all possible lineage sizes.
#'
#' @param x Observed reads.
#' @param b Birth rate.
#' @param d Death rate.
#' @param t Time.
#' @param K Total count.
#' @param nmax Maximum lineage size.
#' @param Nt Total population size.
#' @param variance_vector Precomputed variance vector.
#' @param log Logical, whether to return the log probability.
#' @return Marginal probability of observing reads.
#' @examples
#' marginal_probability_reads(5, 0.1, 0.05, 10, 100, 50, 1000, rep(1, 51), FALSE)
marginal_probability_reads <- function(x, b, d, t, K, nmax, Nt, variance_vector, log = FALSE){
  pns <- sapply(0:nmax, probability_lineage_size, b, d, t)
  pkns <- sapply(0:nmax, probability_lineage_reads, x, b, d, t, K, nmax, Nt, variance_vector)
  prob_k <- sum(pns * pkns)
  if(log){
    return(log(prob_k))
  } else {
    return(prob_k)
  }
}
