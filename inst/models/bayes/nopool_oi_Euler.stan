// Observation error model with pairmeans ODE solution
// yep.

data {
  
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  real err_obs_sigma_min;
  real err_obs_sigma_max;
  
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
  // real DO_mod_1;
  
  real <lower=err_obs_sigma_min, upper=err_obs_sigma_max> err_obs_sigma;
  
}

// manual says: evaluated once per leapfrog step
transformed parameters {
  
  // Declare temporary variables
  vector [n] DO_mod;

  // Model DO time series
  DO_mod[1] <- DO_obs[1]; // DO_mod_1;
  for(i in 1:(n-1)) {
    DO_mod[i+1] <- (
      DO_mod[i] +
        GPP_daily * coef_GPP[i] + 
        ER_daily * coef_ER[i] + 
        K600_daily * coef_K600[i] * (DO_sat_pairmean[i] - DO_mod[i]/2)
    ) / (1 + K600_daily * coef_K600[i] /2);
  }
  
}

model {
  
  
  // Observation error: Compare all the DO predictions to their observations
  DO_obs ~ normal(DO_mod, err_obs_sigma);
  // Prior on the sd of the errors between modeled and observed DO; this is actually unnecessary to specify
  //err_obs_sigma ~ uniform(err_obs_sigma_min, err_obs_sigma_max);
  
  // Daily mean values of GPP and ER (gO2 m^-2 d^-1) and K600 (m d^-1)
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
  
}
