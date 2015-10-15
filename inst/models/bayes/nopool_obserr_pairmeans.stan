// models/bayes/stan/nopool_obserr_pairmeans.txt

data {
  
  # hyperparameters
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  real err_obs_sigma_min;
  real err_obs_sigma_max;
  
  int <lower=0> n; // number of observations in the day
  
  vector [n] DO_obs;
  vector [n] DO_sat;
  vector [n] frac_GPP; // fractions, summing to 1, to partition GPP_daily into per-timestep rates
  vector [n] frac_ER; // fractions, summing to 1, to partition ER_daily into per-timestep rates
  vector [n] frac_D; // fractions, summing to 1, to partition K600_daily into per-timestep rates
  vector [n] depth;
  vector [n] KO2_conv;
  
}

parameters {
  
  real GPP_daily;
  real ER_daily;
  real K600_daily;
  
  // bounded prior (implied uniform) on the sd of the errors between modeled and observed DO
  real <lower=err_obs_sigma_min, upper=err_obs_sigma_max> err_obs_sigma;
  
}

model {
  
  // Declare temporary variables
  vector [n] GPP;
  vector [n] ER;
  vector [n] K;
  vector [n] DO_mod;
  
  // Convert daily rates to per-observation rates
  for(i in 2:n) {
    GPP[i] <- GPP_daily * (frac_GPP[i]+frac_GPP[i-1])/2 / ((depth[i]+depth[i-1])/2);
    ER[i] <- ER_daily * (frac_ER[i]+frac_ER[i-1])/2 / ((depth[i]+depth[i-1])/2);
    K[i] <- K600_daily * (KO2_conv[i]+KO2_conv[i-1])/2 * (frac_D[i]+frac_D[i-1])/2;
  }
  
  // Model DO time series
  DO_mod[1] <- DO_obs[1];
  for(i in 2:n) {
    DO_mod[i] <- (
      DO_mod[i-1] +
        GPP[i] + 
        ER[i] + 
        K[i] * (DO_sat[i] + DO_sat[i-1] - DO_mod[i-1])/2
    ) / (1 + K[i]/2);
  }
  
  // Observation error: Compare all the DO predictions to their observations
  for (i in 2:n) {
    DO_obs[i] ~ normal(DO_mod[i], err_obs_sigma);
  }
  // Prior on the sd of the errors between modeled and observed DO; this is actually unnecessary to specify
  err_obs_sigma ~ uniform(err_obs_sigma_min, err_obs_sigma_max);
  
  // Daily mean values of GPP and ER (gO2 m^-2 d^-1) and K600 (m d^-1)
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
  
}
