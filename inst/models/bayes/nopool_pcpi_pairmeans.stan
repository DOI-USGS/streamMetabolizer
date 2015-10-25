data {
  
  # hyperparameters
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
  real err_proc_iid_sigma_min;
  real err_proc_iid_sigma_max;
  
  int <lower=0> n; // number of observations in the day
  
  vector [n] DO_obs;
  vector [n] DO_sat;
  vector [n] frac_GPP; // fractions, summing to 1, to partition GPP_daily into per-timestep rates
  vector [n] frac_ER; // fractions, summing to 1, to partition ER_daily into per-timestep rates
  vector [n] frac_D; // fractions, summing to 1, to partition K600_daily into per-timestep rates
  vector [n] depth;
  vector [n] KO2_conv;
  
}

// manual says: the statements in the transformed data block are only ever evaluated once
transformed data {
  
  vector [n] coef_GPP;
  vector [n] coef_ER;
  vector [n] coef_K600;
  vector [n] DO_sat_pairmean;
  vector [n] dDOdt_obs;
  vector [n] DO_obs_prev;

  // Prepare vectors of effective length n-1. The first observation is never 
  // actually used but must be specified, so is set to 0
  coef_GPP[1]  <- 0; //frac_GPP[1] / depth[1];
  coef_ER[1]   <- 0; //frac_ER[1] / depth[1];
  coef_K600[1] <- 0; //KO2_conv[1] * frac_D[1];
  DO_sat_pairmean[1] <- 0; //DO_sat[1];
  dDOdt_obs[1] <- 0;
  DO_obs_prev[1] <- 0; //DO_obs[1];
  for(i in 2:n) {
    # Daily to per-timestep conversions using pairmeans
    coef_GPP[i]  <- (frac_GPP[i]+frac_GPP[i-1])/2 / ((depth[i]+depth[i-1])/2);
    coef_ER[i]   <- (frac_ER[i]+frac_ER[i-1])/2 / ((depth[i]+depth[i-1])/2);
    coef_K600[i] <- (KO2_conv[i]+KO2_conv[i-1])/2 * (frac_D[i]+frac_D[i-1])/2;
    DO_sat_pairmean[i] <- (DO_sat[i] + DO_sat[i-1])/2;
    # Not pairmeans, but also needs the 2:n loop: compute stepwise differences in DO
    dDOdt_obs[i] <- (DO_obs[i] - DO_obs[i-1])/frac_D[i]; // assumes frac_D is timestep in days
    # Not pairmeans, but permits a vectorized computation of dDOdt_mod
    DO_obs_prev[i] <- DO_obs[i-1];
  }

}

parameters {
  
  real GPP_daily;
  real ER_daily;
  real K600_daily;
  vector [n] err_proc_acor;
  
  // bounded prior (implied uniform) on the sd of the errors between modeled and observed DO
  real <lower=err_proc_acor_phi_min, upper=err_proc_acor_phi_max> err_proc_acor_phi;
  real <lower=err_proc_acor_sigma_min, upper=err_proc_acor_sigma_max> err_proc_acor_sigma;
  real <lower=err_proc_iid_sigma_min, upper=err_proc_iid_sigma_max> err_proc_iid_sigma;
  
}

// manual says: evaluated once per leapfrog step
transformed parameters {
  
  // Declare temporary variables
  vector [n] GPP;
  vector [n] ER;
  vector [n] K;
  vector [n] dDOdt_mod;
  
  // Convert rates to per-observation rates (vectorized)
  GPP <- GPP_daily * coef_GPP;
  ER <- ER_daily * coef_ER;
  K <- K600_daily * coef_K600;
  
  // Model DO time series; vectorized
  dDOdt_mod <- (
    GPP_daily * coef_GPP + 
    ER_daily * coef_ER + 
    (K600_daily * coef_K600) .* (DO_sat_pairmean - DO_obs_prev/2.0) + 
    err_proc_acor
  ) ./ (1.0 + (K600_daily * coef_K600)/2.0);

}

model {
  
  
  // Autocorrelated process error: Build an error timeseries with a fitted autocorrelation structure
  // err_proc_acor[1] ~ normal(0, err_proc_acor_sigma);
  for(i in 2:n) {
    err_proc_acor[i] ~ normal(err_proc_acor_phi*err_proc_acor[i-1], err_proc_acor_sigma);
  }
  // Prior on the autocorrelation & sd of the process errors
  //err_proc_acor_phi ~ uniform(err_proc_acor_phi_min, err_proc_acor_phi_max); // implied in parameters declaration
  //err_proc_acor_sigma ~ uniform(err_proc_acor_sigma_min, err_proc_acor_sigma_max); // implied in parameters declaration
  
  // Uncorrelated (IID) process error: Compare all the dDOdt predictions to their observations
  dDOdt_obs ~ normal(dDOdt_mod, err_proc_iid_sigma);
  // Prior on the sd of the errors between modeled and observed DO; this is actually unnecessary to specify
  //err_proc_iid_sigma ~ uniform(err_proc_iid_sigma_min, err_proc_iid_sigma_max); // implied in parameters declaration
  
  // Daily mean values of GPP and ER (gO2 m^-2 d^-1) and K600 (m d^-1)
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
  
}