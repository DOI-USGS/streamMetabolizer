// b_Kn_pc_eu_plrcko.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  
  // Parameters of hierarchical priors on K600_daily (normal model)
  real K600_daily_mu_mu;
  real K600_daily_mu_sigma;
  real K600_daily_sigma_location;
  real K600_daily_sigma_scale;
  
  // Error distributions
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
  vector[d] coef_K600_full[n-1];
  vector[d] dDO_obs[n-1];
  
  for(i in 1:(n-1)) {
    // Coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)
    coef_GPP[i]  <- frac_GPP[i] ./ depth[i];
    coef_ER[i]   <- frac_ER[i] ./ depth[i];
    coef_K600_full[i] <- KO2_conv[i] .* frac_D[i] .*
      (DO_sat[i] - DO_obs[i]);
    // dDO observations
    dDO_obs[i] <- DO_obs[i+1] - DO_obs[i];
  }
}

parameters {
  vector[d] GPP_daily;
  vector[d] ER_daily;
  vector<lower=0>[d] K600_daily;
  
  real K600_daily_mu;
  real<lower=0> K600_daily_sigma;
  
  vector[d] err_proc_acor_inc[n-1];
  
  real<lower=0, upper=1> err_proc_acor_phi;
  real<lower=0> err_proc_acor_sigma;
}

transformed parameters {
  vector[d] dDO_mod[n-1];
  vector[d] err_proc_acor[n-1];
  
  // Model DO time series
  // * Euler version
  // * no observation error
  // * autocorrelated process error
  // * reaeration depends on DO_obs
  
  err_proc_acor[1] <- err_proc_acor_inc[1];
  for(i in 1:(n-2)) {
    err_proc_acor[i+1] <- err_proc_acor_phi * err_proc_acor[i] + err_proc_acor_inc[i+1];
  }
  
  // dDO model
  for(i in 1:(n-1)) {
    dDO_mod[i] <- 
      err_proc_acor[i] +
      GPP_daily  .* coef_GPP[i] +
      ER_daily   .* coef_ER[i] +
      K600_daily .* coef_K600_full[i];
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
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);

  // Hierarchical constraints on K600_daily (normal model)
  K600_daily_mu ~ normal(K600_daily_mu_mu, K600_daily_mu_sigma);
  K600_daily_sigma ~ lognormal(K600_daily_sigma_location, K600_daily_sigma_scale);
}
