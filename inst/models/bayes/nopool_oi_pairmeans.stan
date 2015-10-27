// Observation error model with pairmeans ODE solution
// 

data {
  
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  real err_obs_iid_sigma_min;
  real err_obs_iid_sigma_max;
  
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

  // Convert daily rates to per-observation rates using pairmeans. The first 
  // observation is never actually used but must be specified, so is set to 0
  coef_GPP[1]  <- 0; //frac_GPP[1] / depth[1];
  coef_ER[1]   <- 0; //frac_ER[1] / depth[1];
  coef_K600[1] <- 0; //KO2_conv[1] * frac_D[1];
  DO_sat_pairmean[1] <- 0; //DO_sat[1];
  for(i in 2:n) {
    coef_GPP[i]  <- (frac_GPP[i]+frac_GPP[i-1])/2 / ((depth[i]+depth[i-1])/2);
    coef_ER[i]   <- (frac_ER[i]+frac_ER[i-1])/2 / ((depth[i]+depth[i-1])/2);
    coef_K600[i] <- (KO2_conv[i]+KO2_conv[i-1])/2 * (frac_D[i]+frac_D[i-1])/2;
    DO_sat_pairmean[i] <- (DO_sat[i] + DO_sat[i-1])/2;
  }
  
}

parameters {
  
  real GPP_daily;
  real ER_daily;
  real K600_daily;
  // real DO_mod_1;
  
  // bounded prior (implied uniform) on the sd of the errors between modeled and observed DO
  real <lower=err_obs_iid_sigma_min, upper=err_obs_iid_sigma_max> err_obs_iid_sigma;
  
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
    ) / (1 + K[i]/2);
  }

}

model {
  
  
  // Observation error: Compare all the DO predictions to their observations
  DO_obs ~ normal(DO_mod, err_obs_iid_sigma);
  // Prior on the sd of the errors between modeled and observed DO; this is actually unnecessary to specify
  //err_obs_iid_sigma ~ uniform(err_obs_iid_sigma_min, err_obs_iid_sigma_max);
  
  // Daily mean values of GPP and ER (gO2 m^-2 d^-1) and K600 (m d^-1)
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
  
}
