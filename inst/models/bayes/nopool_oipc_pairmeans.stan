// Full state space model for a single day. Matches most closely with
// nopool_pcpi_pairmeans.stan and nopool_oipc_Euler.jags

data {
  
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  real err_proc_acor_phi_min;
  real err_proc_acor_phi_max;
  real err_proc_acor_sigma_min;
  real err_proc_acor_sigma_max;
  real err_obs_iid_sigma_min;
  real err_obs_iid_sigma_max;
  
  int <lower=0> n;
  
  vector [n] DO_obs;
  vector [n] DO_sat;
  vector [n] frac_GPP;
  vector [n] frac_ER;
  vector [n] frac_D;
  vector [n] depth;
  vector [n] KO2_conv;
  
}

// manual says: the statements in the transformed data block are only ever evaluated once
transformed data {
  
  vector [n-1] coef_GPP;
  vector [n-1] coef_ER;
  vector [n-1] coef_K600;
  vector [n-1] DO_sat_pairmean;

  // Convert daily rates to per-observation rates using pairmeans. The first 
  // observation is never actually used but must be specified, so is set to 0
//   for(i in 2:n) {
//     coef_GPP[i]  <- (frac_GPP[i]+frac_GPP[i-1])/2 / ((depth[i]+depth[i-1])/2);
//     coef_ER[i]   <- (frac_ER[i]+frac_ER[i-1])/2 / ((depth[i]+depth[i-1])/2);
//     coef_K600[i] <- (KO2_conv[i]+KO2_conv[i-1])/2 * (frac_D[i]+frac_D[i-1])/2;
//     DO_sat_pairmean[i] <- (DO_sat[i] + DO_sat[i-1])/2;
//   }
  
  // coefficients by pairmeans (e.g., (frac_GPP[i] + frac_GPP[i+1])/2 applies to the DO step from i to i+1)
  for(i in 1:(n-1)) {
    coef_GPP[i]  <- (frac_GPP[i] + frac_GPP[i+1])/2 / ((depth[i] + depth[i+1])/2);
    coef_ER[i]   <- (frac_ER[ i] + frac_ER[ i+1])/2 / ((depth[i] + depth[i+1])/2);
    coef_K600[i] <- (KO2_conv[i] + KO2_conv[i+1])/2 * (frac_D[i] + frac_D[i+1])/2;
    DO_sat_pairmean[i] <- (DO_sat[i] + DO_sat[i+1])/2;
  }
  
}

parameters {
  
  real GPP_daily;
  real ER_daily;
  real K600_daily;
  real <lower=err_proc_acor_phi_min,   upper=err_proc_acor_phi_max>   err_proc_acor_phi;
  real <lower=err_proc_acor_sigma_min, upper=err_proc_acor_sigma_max> err_proc_acor_sigma;
  real <lower=err_proc_iid_sigma_min,  upper=err_proc_iid_sigma_max>  err_proc_iid_sigma;
  vector[n-1] err_proc_acor;
  // real DO_mod_1;
  
}

// manual says: evaluated once per leapfrog step
transformed parameters {
  
  vector [n] DO_mod;
  
  // Model DO time series
  DO_mod[1] <- DO_obs[1]; // DO_mod_1;
  for(i in 2:n) {
    DO_mod[i] <- (
      DO_mod[i-1] +
        GPP_daily * coef_GPP[i] + 
        ER_daily * coef_ER[i] + 
        (K600_daily * coef_K600[i]) .* (DO_sat_pairmean[i] - DO_mod[i-1]/2)
    ) / (1 + K600_daily * coef_K600[i] / 2);
  }

}

model {
  
  // Process error: Build an error timeseries with a fitted autocorrelation structure
  // err_proc_acor[1] ~ normal(0, err_proc_acor_sigma);
  for(i in 2:n) {
    err_proc_acor[i] ~ normal(err_proc_acor_phi*err_proc_acor[i-1], err_proc_acor_sigma);
  }
  // Prior on the autocorrelation & sd of the process errors
  //err_proc_acor_phi ~ uniform(err_proc_acor_phi_min, err_proc_acor_phi_max); // implied in parameters declaration
  //err_proc_acor_sigma ~ uniform(err_proc_acor_sigma_min, err_proc_acor_sigma_max); // implied in parameters declaration
  
  // Observation error: Compare all the DO predictions to their observations
  //DO_mod_1 ~ normal(DO_obs[1], err_obs_iid_sigma);
  DO_obs ~ normal(DO_mod, err_obs_iid_sigma);
  // Prior on the sd of the errors between modeled and observed DO; this is actually unnecessary to specify
  //err_obs_iid_sigma ~ uniform(err_obs_iid_sigma_min, err_obs_iid_sigma_max); // implied in parameters declaration
  
  // Daily mean values of GPP and ER (gO2 m^-2 d^-1) and K600 (m d^-1)
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
  
}
