// b_np_oipc_pm_plrckm.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  // Error distributions
  real err_obs_iid_sigma_location;
  real err_obs_iid_sigma_scale;
  real err_proc_acor_phi_alpha;
  real err_proc_acor_phi_beta;
  real err_proc_acor_sigma_location;
  real err_proc_acor_sigma_scale;
  
  // Data dimensions
  int<lower=1> d; # number of dates
  int<lower=1> n; # number of observations per date
  
  // Daily data
  vector[d] DO_obs_1;
  
  // Data
  vector[d] DO_obs[n];
  vector[d] DO_sat[n];
  vector[d] frac_GPP[n];
  vector[d] frac_ER[n];
  vector[d] frac_D[n];
  vector[d] depth[n];
  vector[d] KO2_conv[n];
}

transformed data {
  vector[d] coef_GPP[n-1];
  vector[d] coef_ER[n-1];
  vector[d] coef_K600_part[n-1];
  vector[d] DO_sat_pairmean[n-1];
  
  for(i in 1:(n-1)) {
    // Coefficients by pairmeans (e.g., mean(frac_GPP[i:(i+1)]) applies to the DO step from i to i+1)
    coef_GPP[i]  <- (frac_GPP[i] + frac_GPP[i+1])/2.0 ./ ((depth[i] + depth[i+1])/2.0);
    coef_ER[i]   <- (frac_ER[i] + frac_ER[i+1])/2.0 ./ ((depth[i] + depth[i+1])/2.0);
    coef_K600_part[i] <- (KO2_conv[i] + KO2_conv[i+1])/2.0 .* (frac_D[i] + frac_D[i+1])/2.0;
    DO_sat_pairmean[i] <- (DO_sat[i] + DO_sat[i+1])/2.0;
  }
}

parameters {
  vector[d] GPP_daily;
  vector[d] ER_daily;
  vector<lower=0>[d] K600_daily;
  
  vector[d] err_proc_acor_inc[n-1];
  
  real<lower=0> err_obs_iid_sigma;
  real<lower=0, upper=1> err_proc_acor_phi;
  real<lower=0> err_proc_acor_sigma;
}

transformed parameters {
  vector[d] DO_mod[n];
  vector[d] err_proc_acor[n-1];
  
  // Model DO time series
  // * pairmeans version
  // * observation error
  // * autocorrelated process error
  // * reaeration depends on DO_mod
  
  err_proc_acor[1] <- err_proc_acor_inc[1];
  for(i in 1:(n-2)) {
    err_proc_acor[i+1] <- err_proc_acor_phi * err_proc_acor[i] + err_proc_acor_inc[i+1];
  }
  
  // DO model
  DO_mod[1] <- DO_obs_1;
  for(i in 1:(n-1)) {
    DO_mod[i+1] <- (
      DO_mod[i] +
      err_proc_acor[i] +
      GPP_daily .* coef_GPP[i] +
      ER_daily .* coef_ER[i] +
      K600_daily .* coef_K600_part[i] .* (DO_sat_pairmean[i] - DO_mod[i]/2.0)
    ) ./ (1.0 + K600_daily .* coef_K600_part[i] / 2.0);
  }
}

model {
  // Autocorrelated process error
  for(i in 1:(n-1)) {
    err_proc_acor_inc[i] ~ normal(0, err_proc_acor_sigma);
  }
  // Autocorrelation (phi) & SD (sigma) of the process errors
  err_proc_acor_phi ~ beta(err_proc_acor_phi_alpha, err_proc_acor_phi_beta);
  err_proc_acor_sigma ~ lognormal(err_proc_acor_sigma_location, err_proc_acor_sigma_scale);
  
  // Independent, identically distributed observation error
  for(i in 1:n) {
    DO_obs[i] ~ normal(DO_mod[i], err_obs_iid_sigma);
  }
  // SD (sigma) of the observation errors
  err_obs_iid_sigma ~ lognormal(err_obs_iid_sigma_location, err_obs_iid_sigma_scale);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
}
