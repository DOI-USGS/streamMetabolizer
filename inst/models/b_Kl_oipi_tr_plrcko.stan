// b_Kl_oipi_tr_plrcko.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  
  // Parameters of hierarchical priors on K600_daily (linear model)
  vector[2] K600_daily_beta_mu;
  vector[2] K600_daily_beta_sigma;
  real<lower=0> K600_daily_sigma_scale;
  
  // Error distributions
  real<lower=0> err_obs_iid_sigma_scale;
  real<lower=0> err_proc_iid_sigma_scale;
  
  // Data dimensions
  int<lower=1> d; # number of dates
  int<lower=1> n; # number of observations per date
  
  // Daily data
  vector[d] DO_obs_1;
  vector[d] ln_discharge_daily;
  
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
    // Coefficients for trapezoid rule (e.g., mean(frac_GPP[i:(i+1)]) applies to the DO step from i to i+1)
    coef_GPP[i] = ((frac_GPP[i] ./ depth[i]) + (frac_GPP[i+1] ./ depth[i+1]))/2.0;
    coef_ER[i] = ((frac_ER[i] ./ depth[i]) + (frac_ER[i+1] ./ depth[i+1]))/2.0;
    coef_K600_full[i] = ((KO2_conv[i] .* frac_D[i]) + (KO2_conv[i+1] .* frac_D[i+1]))/2.0 .*
      (DO_sat[i] + DO_sat[i+1] - DO_obs[i] - DO_obs[i+1])/2.0;
    // dDO observations
    dDO_obs[i] = DO_obs[i+1] - DO_obs[i];
  }
}

parameters {
  vector[d] GPP_daily;
  vector[d] ER_daily;
  vector<lower=0>[d] K600_daily;
  
  vector[2] K600_daily_beta;
  real<lower=0> K600_daily_sigma_scaled;
  
  real<lower=0> err_obs_iid_sigma_scaled;
  real<lower=0> err_proc_iid_sigma_scaled;
}

transformed parameters {
  real<lower=0> K600_daily_sigma;
  vector[d] K600_daily_pred;
  real<lower=0> err_obs_iid_sigma;
  real<lower=0> err_proc_iid_sigma;
  vector[d] DO_mod[n];
  vector[d] dDO_mod[n-1];
  
  // Rescale pooling & error distribution parameters
  K600_daily_sigma = K600_daily_sigma_scale * K600_daily_sigma_scaled;
  err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled;
  err_proc_iid_sigma = err_proc_iid_sigma_scale * err_proc_iid_sigma_scaled;
  
  // Hierarchical, linear model of K600_daily
  K600_daily_pred = K600_daily_beta[1] + K600_daily_beta[2] * ln_discharge_daily;
  
  // Model DO time series
  // * trapezoid version
  // * observation error
  // * IID process error
  // * reaeration depends on DO_obs
  
  // dDO model
  for(i in 1:(n-1)) {
    dDO_mod[i] = 
      GPP_daily  .* coef_GPP[i] +
      ER_daily   .* coef_ER[i] +
      K600_daily .* coef_K600_full[i];
  }
  
  // DO model
  DO_mod[1] = DO_obs_1;
  for(i in 1:(n-1)) {
    DO_mod[i+1] = (
      DO_mod[i] +
      dDO_mod[i] + err_proc_iid[i]);
  }
}

model {
  // Process error
  for(i in 1:(n-1)) {
    // Independent, identically distributed process error
    dDO_obs[i] ~ normal(dDO_mod[i], err_proc_iid_sigma);
  }
  // SD (sigma) of the IID process errors
  err_proc_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Independent, identically distributed observation error
  for(i in 2:n) {
    DO_obs[i] ~ normal(DO_mod[i], err_obs_iid_sigma);
  }
  // SD (sigma) of the observation errors
  err_obs_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_pred, K600_daily_sigma);

  // Hierarchical constraints on K600_daily (linear model)
  K600_daily_beta ~ normal(K600_daily_beta_mu, K600_daily_beta_sigma);
  K600_daily_sigma_scaled ~ cauchy(0, 1);
}
