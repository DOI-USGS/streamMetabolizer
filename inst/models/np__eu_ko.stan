// np__eu_ko.stan

data {
  // Metabolism distributions
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  // Error distributions
  
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
  vector [n-1] dDO_obs;
  
  for(i in 1:(n-1)) {
    // Coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)
    coef_GPP[i]  <- frac_GPP[i] / depth[i];
    coef_ER[i]   <- frac_ER[ i] / depth[i];
    coef_K600_full[i] <- KO2_conv[i] * frac_D[i] * 
      (DO_sat[i] - DO_obs[i]);
    // dDO observations
    dDO_obs[i] <- DO_obs[i+1] - DO_obs[i];
  }
}

parameters {
  real GPP_daily;
  real ER_daily;
  real K600_daily;
}

transformed parameters {
  vector [n-1] dDO_mod;
  
  // Model DO time series
  // * Euler version
  // * no observation error
  // * no process error
  // * reaeration depends on DO_obs
  
  // dDO model
  dDO_mod <- 
    GPP_daily * coef_GPP +
    ER_daily * coef_ER +
    K600_daily * coef_K600_full;
}

model {
  // Daily metabolism values
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
}
