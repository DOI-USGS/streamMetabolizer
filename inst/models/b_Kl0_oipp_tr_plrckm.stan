// b_Kl0_oipp_tr_plrckm.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_lower;
  real<lower=0> GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_upper;
  real<lower=0> ER_daily_sigma;
  
  // Parameters of hierarchical priors on K600_daily (linear_sdzero model)
  real lnK600_lnQ_intercept_mu;
  real<lower=0> lnK600_lnQ_intercept_sigma;
  real lnK600_lnQ_slope_mu;
  real<lower=0> lnK600_lnQ_slope_sigma;
  
  // Error distributions
  real<lower=0> err_obs_iid_sigma_scale;
  real<lower=0> err_mult_GPP_sdlog_sigma;
  
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
  
  real lnK600_lnQ_intercept;
  real lnK600_lnQ_slope;
  
  real<lower=0> err_obs_iid_sigma_scaled;
  real<lower=0> err_mult_GPP_sdlog_scaled;
  vector<lower=0>[d] err_mult_GPP[n];
}

transformed parameters {
  vector[d] K600_daily_predlog;
  vector[d] K600_daily;
  real<lower=0> err_obs_iid_sigma;
  real<lower=0> err_mult_GPP_sdlog;
  vector[d] GPP_inst[n];
  vector[d] ER_inst[n];
  vector[d] KO2_inst[n];
  vector<lower=0>[d] combo_mult_GPP[n];
  vector<lower=0>[d] mean_combo_mult_GPP;
  vector[d] DO_mod[n];
  
  // Rescale error distribution parameters
  err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled;
  err_mult_GPP_sdlog = err_mult_GPP_sdlog_sigma * err_mult_GPP_sdlog_scaled;
  
  // Hierarchical, linear_sdzero model of K600_daily
  K600_daily_predlog = lnK600_lnQ_intercept + lnK600_lnQ_slope * lnQ_daily;
  K600_daily = exp(K600_daily_predlog);
  
  // Model DO time series
  // * trapezoid version
  // * observation error
  // * no process error
  // * reaeration depends on DO_mod
  
  // Calculate individual process rates
  for(i in 1:n) {
    combo_mult_GPP[i] = err_mult_GPP[i] .* light_mult_GPP[i];
  }
  for(j in 1:d) {
    mean_combo_mult_GPP[j] = sum(combo_mult_GPP[,j]) / n;
  }
  for(i in 1:n) {
    GPP_inst[i] = GPP_daily .* combo_mult_GPP[i] ./ mean_combo_mult_GPP;
    ER_inst[i] = ER_daily .* const_mult_ER[i];
    KO2_inst[i] = K600_daily .* KO2_conv[i];
  }
  
  // DO model
  DO_mod[1] = DO_obs_1;
  for(i in 1:(n-1)) {
    DO_mod[i+1] =
      DO_mod[i] .*
        (2.0 - KO2_inst[i] * timestep) ./ (2.0 + KO2_inst[i+1] * timestep) + (
        (GPP_inst[i] + ER_inst[i]) ./ depth[i] +
        (GPP_inst[i+1] + ER_inst[i+1]) ./ depth[i+1] +
        KO2_inst[i] .* DO_sat[i] + KO2_inst[i+1] .* DO_sat[i+1]
      ) .* (timestep ./ (2.0 + KO2_inst[i+1] * timestep));
  }
}

model {
  // GPP-only independent, identically distributed process error
  for(i in 1:n) {
    err_mult_GPP[i] ~ lognormal(0, err_mult_GPP_sdlog);
  }
  // SD (sigma) of the GPP IID process error multipliers
  err_mult_GPP_sdlog_scaled ~ normal(0, 1);
  
  // Independent, identically distributed observation error
  for(i in 2:n) {
    DO_obs[i] ~ normal(DO_mod[i], err_obs_iid_sigma);
  }
  // SD (sigma) of the observation errors
  err_obs_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  // Hierarchical constraints on K600_daily (linear_sdzero model)
  lnK600_lnQ_intercept ~ normal(lnK600_lnQ_intercept_mu, lnK600_lnQ_intercept_sigma);
  lnK600_lnQ_slope ~ normal(lnK600_lnQ_slope_mu, lnK600_lnQ_slope_sigma);
  
}
generated quantities {
  vector[d] err_obs_iid[n];
  vector[d] GPP_inst_partial[n];
  vector[d] err_proc_GPP[n];
  int n_light_day; // temporary
  vector[n] GPP_inst_day; // temporary
  vector[n] GPP_inst_diff_day; // temporary
  vector[d] GPP_pseudo_R2;
  vector[d] GPP;
  vector[d] ER;
  vector[n] DO_obs_vec; // temporary
  vector[n] DO_mod_vec; // temporary
  vector[d] DO_R2;
  
  for(i in 1:n) {
    err_obs_iid[i] = DO_mod[i] - DO_obs[i];
  }
  for(i in 1:n) {
    GPP_inst_partial[i] = GPP_daily .* light_mult_GPP[i];
    err_proc_GPP[i] = GPP_inst[i] - GPP_inst_partial[i];
  }
  GPP_inst_day = rep_vector(0, n);
  GPP_inst_diff_day = rep_vector(0, n);
  for(j in 1:d) {
    GPP[j] = sum(GPP_inst[1:n24,j]) / n24;
    ER[j] = sum(ER_inst[1:n24,j]) / n24;
    
    // Compute R2 for DO observations relative to the modeled, process-error-corrected state (DO_mod)
    for(i in 1:n) {
      DO_mod_vec[i] = DO_mod[i,j];
      DO_obs_vec[i] = DO_obs[i,j];
    }
    DO_R2[j] = 1 - sum((DO_mod_vec - DO_obs_vec) .* (DO_mod_vec - DO_obs_vec)) / sum((DO_obs_vec - mean(DO_obs_vec)) .* (DO_obs_vec - mean(DO_obs_vec)));
    
    // Compute GPP_pseudo_R2 (because model has GPP process error)
    n_light_day = 0;
    for(i in 1:n) {
      if(light_mult_GPP[i,j] > 0) {
        n_light_day += 1;
        GPP_inst_day[n_light_day] = GPP_inst[i,j];
        GPP_inst_diff_day[n_light_day] = GPP_inst[i,j] - GPP_inst_partial[i,j];
      }
    }
    GPP_pseudo_R2[j] = 1 - variance(GPP_inst_diff_day[1:n_light_day]) / variance(GPP_inst_day[1:n_light_day]);
  }
  
}
