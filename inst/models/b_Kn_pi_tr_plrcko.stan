// b_Kn_pi_tr_plrcko.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  
  // Parameters of hierarchical priors on K600_daily (normal model)
  real K600_daily_meanlog_meanlog;
  real K600_daily_meanlog_sdlog;
  real<lower=0> K600_daily_sdlog_scale;
  
  // Error distributions
  real<lower=0> err_proc_iid_sigma_scale;
  
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
  real<lower=0> timestep; # length of each timestep in days
  timestep = frac_D[1,1];
}

parameters {
  vector[d] GPP_daily;
  vector[d] ER_daily;
  vector[d] K600_daily;
  
  real K600_daily_mu;
  real<lower=0> K600_daily_sdlog_scaled;
  
  real<lower=0> err_proc_iid_sigma_scaled;
  
  vector[d] err_proc_iid[n-1];
}

transformed parameters {
  real<lower=0> K600_daily_sdlog;
  vector[d] DO_mod_partial_sigma[n];
  real<lower=0> err_proc_iid_sigma;
  vector[d] GPP[n];
  vector[d] ER[n];
  vector[d] KO2[n];
  vector[d] DO_mod_partial[n];
  
  // Rescale pooling & error distribution parameters
  K600_daily_sdlog = K600_daily_sdlog_scale * K600_daily_sdlog_scaled;
  err_proc_iid_sigma = err_proc_iid_sigma_scale * err_proc_iid_sigma_scaled;
  
  // Model DO time series
  // * trapezoid version
  // * no observation error
  // * IID process error
  // * reaeration depends on DO_obs
  
  // Calculate individual process rates
  for(i in 1:n) {
    GPP[i] = GPP_daily .* frac_GPP[i];
    ER[i] = ER_daily .* frac_ER[i];
    KO2[i] = K600_daily .* KO2_conv[i];
  }
  
  // DO model
  for(i in 1:(n-1)) {
    DO_mod_partial[i+1] =
      DO_obs[i] + (
        - KO2[i] .* DO_obs[i] - KO2[i+1] .* DO_obs[i+1] +
        (GPP[i] + ER[i]) ./ depth[i] +
        (GPP[i+1] + ER[i+1]) ./ depth[i+1] +
        KO2[i] .* DO_sat[i] + KO2[i+1] .* DO_sat[i+1]
      ) * (timestep / 2.0);
    for(j in 1:d) {
      DO_mod_partial_sigma[i+1,j] = err_proc_iid_sigma * 
        sqrt(pow(depth[i,j], -2) + pow(depth[i+1,j], -2)) .*
        (timestep / 2.0);
    }
  }
}

model {
  // Process error
  for(i in 2:n) {
    // Independent, identically distributed process error
    DO_obs[i] ~ normal(DO_mod_partial[i], DO_mod_partial_sigma[i]);
  }
  // SD (sigma) of the IID process errors
  err_proc_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ lognormal(K600_daily_pred, K600_daily_sdlog);

  // Hierarchical constraints on K600_daily (normal model)
  K600_daily_pred ~ lognormal(K600_daily_meanlog_meanlog, K600_daily_meanlog_sdlog );
  K600_daily_sdlog_scaled ~ cauchy(0, 1);
}
