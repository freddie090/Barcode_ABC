
functions {
    
    // This function returns the log probability of seeing a barcode n times 
    // given b, d, and t.

    real birthdeath_prob_n_nl(real n, real b, real d, real t){
        real log_pn;
        real alpha;
        real beta;
        
        alpha = d * (exp((b - d) * t) - 1) / (b * exp((b - d) * t) - d);
        beta = b * (exp((b - d) * t) - 1) / (b * exp((b - d) * t) - d);
        if(n == 0){
            log_pn = log(alpha);
        } else {
            log_pn = log(1 - alpha) + log(1 - beta) + (n - 1)*log(beta);
        }
        
    return log_pn;

    }

    // This version conditions on 'survival' - for negative binomial 
    // (over-dispersed Poisson) sampling, that a lineage is sampled > 0 times.
    
    real samp_K_bd_cs_like_lpmf(int x, real b, real d, real t, int K, int nmax, 
                                real Nt, vector jvar_vec){
        
        real lprob;
        vector[nmax + 1] log_pxnns; // log_p(xi | n) + log_p(n)
        vector[nmax + 1] log_pgt0nns; //
        vector[nmax + 1] log_pn; // log_p(n)
        real log_pxn; // log_p(xi | n)
        real ngbin_mu; // mean of negative binomial given n and Nt.
        real ngbin_phi; // precision of negative binomial given var(j').
    
        // Pre-calculate the log(p(n))s to speed up.   
        for(n in 0:nmax){
            // log(p(n ; b, d, t)). Adjust to n + 1 for indexing.
            log_pn[n + 1] = birthdeath_prob_n_nl(n, b, d, t);
        }
        
        // First calculate log_p(k = xi) = log_p(k = xi | n) + log_p(n)
        // summed over all ns. 
        for(n in 0:nmax){

            if(n == 0){
                if(x == 0){
                    log_pxn = log(1);
                } else {
                    log_pxn = log(0);
                }
            } else {
                // mean and precision of negative binomial. 
                ngbin_mu = K*(n/Nt);
                // Convert desired variance into a precision value. 
                ngbin_phi = (ngbin_mu^2)/(jvar_vec[n + 1] - ngbin_mu);
                // log_p(k = xi | n)
                log_pxn = neg_binomial_2_lpmf(x | ngbin_mu, ngbin_phi);
            }
            // Adjust to n + 1 for indexing.
            log_pxnns[n + 1] = log_pn[n + 1] + log_pxn;
        }

    lprob = log_sum_exp(log_pxnns);
        
    return lprob;
        
    }

}

data {

    int<lower=0> N0;       // Number of uniquely barcoded cells at t=0.
    int<lower=0> Nc;       // size of count/frequency vector. 
    int<lower=0> Nvec[Nc]; // 1:size of count/frequency vector. 
    real t;
    int<lower=0> K;
    int<lower=0> I;       // Number of replicates. 
    int<lower=0> nmax;
    real Nt;
    int<lower=0> Ncv[I];  // size of each replicate's count/freq vec.
    int<lower=0> csv[Nc]; // vector of combined counts. 
    int<lower=0> fsv[Nc]; // vector of combined frequencies of the counts.
    matrix<lower=0> [(nmax+1), I] jvar_mat; // vector of var(j'). 

}

transformed data {
}

parameters {
    
    real<lower=0, upper=1> s_scaled;    
    real<lower=0> l;     // s <= l

}

transformed parameters {
    
    real s;
    real b; 
    real d;
    real Nm;
    real Nvar;
    real Nsd;
    
    s = l*s_scaled;

    // equivalent to transform b=(l+s)/2, d=(l-s)/2
    b = (l+s)/2;
    d = (l-s)/2;
    
    // Mean population size given s.
    Nm = N0*(exp(s*t));
    // Std.dev population size given s.
    Nvar = N0*(l/s)*exp(s*t)*(exp(s*t) - 1);
    Nsd = sqrt(Nvar);
    
}

model {

    // To update the replicate location in count/frequency vectors.
    int pos;

    // priors
    s_scaled ~ beta(1, 2);    
    l ~ gamma(1, 0.5);

    // fits
    Nt ~ normal(Nm, Nsd); 
    
    pos = 1; 

    // Compressed counts version. 
    for(i in 1:I){
        for(x in segment(Nvec, pos, Ncv[i])){
            target += fsv[x] * ( samp_K_bd_cs_like_lpmf(csv[x] | b, d, t, K, nmax, Nt, jvar_mat[, i]));
        }
        pos = pos + Ncv[i];
    }

    
}
    
generated quantities {
    
    int ll_pos;
    int Nc_pos;
    vector[Nc] log_lik;
    
    // Vector of log-likelihoods for the given posterior draw values of b and d. 
    // For multiple replicates - just save for every unique count. Don't worry 
    // about segmenting until post-processing in R. 
    ll_pos = 1; 
    Nc_pos = 1;
    
    
    for(i in 1:I){
        for(x in segment(Nvec, ll_pos, Ncv[i])){
            log_lik[Nc_pos] = samp_K_bd_cs_like_lpmf(csv[x] | b, d, t, K, nmax, Nt, jvar_mat[, i]);
            Nc_pos = Nc_pos + 1;
        }
        ll_pos = ll_pos + Ncv[i];
    }
    

}
