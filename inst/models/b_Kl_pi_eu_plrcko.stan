// b_Kl_pi_eu_plrcko.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_lower;
  real<lower=0> GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_upper;
  real<lower=0> ER_daily_sigma;
  
  // Parameters of hierarchical priors on K600_daily (linear model)
  real lnK600_lnQ_intercept_mu;
  real<lower=0> lnK600_lnQ_intercept_sigma;
  real lnK600_lnQ_slope_mu;
  real<lower=0> lnK600_lnQ_slope_sigma;
  real<lower=0> K600_daily_sigma_sigma;
  
  // Error distributions
  real<lower=0> err_proc_iid_sigma_scale;
  
  // Data dimensions
  int<lower=1> d; // number of dates
  real<lower=0> timestep; // length of each timestep in days
  int<lower=1> n24; // number of observations in first 24 hours per date
  int<lower=1> n; // number of observations per date
  
  // Daily data
  vector[d] DO_obs_1;
  vector[d] lnQ_daily;
  
  // Data
  vector[d] DO_obs[n];
  vector[d] DO_sat[n];
  vector[d] light_mult_GPP[n];
  vector[d] const_mult_ER[n];
  vector[d] depth[n];
  vector[d] KO2_conv[n];
}

parameters {
  vector<lower=GPP_daily_lower>[d] GPP_daily;
  vector<upper=ER_daily_upper>[d] ER_daily;
  vector<lower=0>[d] K600_daily;
  
  real lnK600_lnQ_intercept;
  real lnK600_lnQ_slope;
  real<lower=0> K600_daily_sigma_scaled;
  
  real<lower=0> err_proc_iid_sigma_scaled;
}

transformed parameters {
  vector[d] K600_daily_predlog;
  real<lower=0> K600_daily_sigma;
  vector[d] DO_mod_partial_sigma[n];
  real<lower=0> err_proc_iid_sigma;
  vector[d] GPP_inst[n];
  vector[d] ER_inst[n];
  vector[d] KO2_inst[n];
  vector[d] DO_mod_partial[n];
  
  // Rescale pooling distribution parameter
  K600_daily_sigma = K600_daily_sigma_sigma * K600_daily_sigma_scaled;
  
  // Rescale error distribution parameters
  err_proc_iid_sigma = err_proc_iid_sigma_scale * err_proc_iid_sigma_scaled;
  
  // Hierarchical, linear model of K600_daily
  K600_daily_predlog = lnK600_lnQ_intercept + lnK600_lnQ_slope * lnQ_daily;
  
  // Model DO time series
  // * euler version
  // * no observation error
  // * IID process error
  // * reaeration depends on DO_obs
  
  // Calculate individual process rates
  for(i in 1:n) {
    GPP_inst[i] = GPP_daily .* light_mult_GPP[i];
    ER_inst[i] = ER_daily .* const_mult_ER[i];
    KO2_inst[i] = K600_daily .* KO2_conv[i];
  }
  
  // DO model
  DO_mod_partial[1] = DO_obs_1;
  DO_mod_partial_sigma[1] = err_proc_iid_sigma * timestep ./ depth[1];
  for(i in 1:(n-1)) {
    DO_mod_partial[i+1] =
      DO_obs[i] + (
        (GPP_inst[i] + ER_inst[i]) ./ depth[i] +
        KO2_inst[i] .* (DO_sat[i] - DO_obs[i])
      ) * timestep;
    for(j in 1:d) {
      DO_mod_partial_sigma[i+1,j] = err_proc_iid_sigma * 
        timestep ./ depth[i,j];
    }
  }
}

model {
  // Independent, identically distributed process error
  for(i in 1:n) {
    DO_obs[i] ~ normal(DO_mod_partial[i], DO_mod_partial_sigma[i]);
  }
  // SD (sigma) of the IID process errors
  err_proc_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(exp(K600_daily_predlog), K600_daily_sigma);
  // Hierarchical constraints on K600_daily (linear model)
  lnK600_lnQ_intercept ~ normal(lnK600_lnQ_intercept_mu, lnK600_lnQ_intercept_sigma);
  lnK600_lnQ_slope ~ normal(lnK600_lnQ_slope_mu, lnK600_lnQ_slope_sigma);
  K600_daily_sigma_scaled ~ normal(0, 1);
  
}
generated quantities {
  vector[d] err_proc_iid[n-1];
  vector[d] GPP;
  vector[d] ER;
  vector[d] DO_R2;
  
  for(i in 2:n) {
    err_proc_iid[i-1] = (DO_mod_partial[i] - DO_obs[i]) .* (err_proc_iid_sigma ./ DO_mod_partial_sigma[i]);
  }
  for(j in 1:d) {
    GPP[j] = sum(GPP_inst[1:n24,j]) / n24;
    ER[j] = sum(ER_inst[1:n24,j]) / n24;
    
    // R2 for DO observations is always 1 for process-error-only models
    DO_R2[j] = 1;
  }
  
}
