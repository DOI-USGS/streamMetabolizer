// np_oi_eu_km.stan

data {

  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
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

transformed data {

  vector [n-1] coef_GPP;
  vector [n-1] coef_ER;
  vector [n-1] coef_K600_part;
  vector [n] DO_obs_copy;
  
  for(i in 1:n) {
    DO_obs_copy[i] <- DO_obs[i];
  }
  
  // Coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)
  for(i in 1:(n-1)) {
    coef_GPP[i]  <- frac_GPP[i] / depth[i];
    coef_ER[i]   <- frac_ER[ i] / depth[i];
    coef_K600_part[i] <- KO2_conv[i] * frac_D[i];
  }
  
}

parameters {

  real GPP_daily;
  real ER_daily;
  real K600_daily;
  
  real <lower=err_obs_iid_sigma_min, upper=err_obs_iid_sigma_max> err_obs_iid_sigma;
  
}

transformed parameters {

  vector [n] DO_mod;
  
  // Model DO time series (Euler version)
  DO_mod[1] <- DO_obs_copy[1];
  for(i in 1:(n-1)) {
    DO_mod[i+1] <- (
      DO_mod[i] +
      GPP_daily * coef_GPP[i] +
      ER_daily * coef_ER[i] +
      K600_daily * coef_K600_part[i] * (DO_sat[i] - DO_mod[i])
    );
  }
  
}

model {

  // Priors on observation error
  for(i in 1:n) {
    DO_obs[i] ~ normal(DO_mod[i], err_obs_iid_sigma);
  }
  err_obs_iid_sigma ~ uniform(err_obs_iid_sigma_min, err_obs_iid_sigma_max);
  
  // Priors on daily metabolism values
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
  
}
