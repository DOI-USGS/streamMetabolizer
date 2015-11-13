// b_np_oipi_pm_ko.stan

data {
  // Metabolism distributions
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  // Error distributions
  real err_obs_iid_sigma_min;
  real err_obs_iid_sigma_max;
  real err_proc_iid_sigma_min;
  real err_proc_iid_sigma_max;
  
  // Daily data
  int <lower=0> n;
  real DO_obs_1;
  
  // Data
  vector [n] DO_obs;
  vector [n] DO_sat;
  vector [n] frac_GPP;
  vector [n] frac_ER;
  vector [n] frac_D;
  vector [n] depth;
  vector [n] KO2_conv;
}

transformed data {
  vector [n-1] coef_GPP;
  vector [n-1] coef_ER;
  vector [n-1] coef_K600_full;
  
  for(i in 1:(n-1)) {
    // Coefficients by pairmeans (e.g., mean(frac_GPP[i:(i+1)]) applies to the DO step from i to i+1)
    coef_GPP[i]  <- (frac_GPP[i] + frac_GPP[i+1])/2 / ((depth[i] + depth[i+1])/2);
    coef_ER[i]   <- (frac_ER[ i] + frac_ER[ i+1])/2 / ((depth[i] + depth[i+1])/2);
    coef_K600_full[i] <- (KO2_conv[i] + KO2_conv[i+1])/2 * (frac_D[i] + frac_D[i+1])/2 *
      (DO_sat[i] + DO_sat[i+1] - DO_obs[i] - DO_obs[i+1])/2;
  }
}

parameters {
  real GPP_daily;
  real ER_daily;
  real K600_daily;
  
  real <lower=err_obs_iid_sigma_min,   upper=err_obs_iid_sigma_max>  err_obs_iid_sigma;
  real <lower=err_proc_iid_sigma_min,  upper=err_proc_iid_sigma_max>  err_proc_iid_sigma;
}

transformed parameters {
  vector [n] DO_mod;
  vector [n-1] dDO_mod;
  
  // Model DO time series
  // * pairmeans version
  // * observation error
  // * IID process error
  // * reaeration depends on DO_obs
  
  // dDO model
  dDO_mod <- 
    GPP_daily * coef_GPP +
    ER_daily * coef_ER +
    K600_daily * coef_K600_full;
  
  // DO model
  DO_mod[1] <- DO_obs_1;
  for(i in 1:(n-1)) {
    DO_mod[i+1] <- (
      DO_mod[i] +
      dDO_mod[i]);
  }
}

model {
  // Independent, identically distributed process error
  for (i in 1:(n-1)) {
    dDO_obs[i] ~ normal(dDO_mod[i], err_proc_iid_sigma);
  }
  err_proc_iid_sigma ~ uniform(err_proc_iid_sigma_min, err_proc_iid_sigma_max);
  
  // Independent, identically distributed observation error
  for(i in 1:n) {
    DO_obs[i] ~ normal(DO_mod[i], err_obs_iid_sigma);
  }
  // SD (sigma) of the observation errors
  err_obs_iid_sigma ~ uniform(err_obs_iid_sigma_min, err_obs_iid_sigma_max);
  
  // Daily metabolism values
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
}
