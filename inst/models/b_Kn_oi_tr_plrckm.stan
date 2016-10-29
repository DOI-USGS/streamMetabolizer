// b_Kn_oi_tr_plrckm.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  
  // Parameters of hierarchical priors on K600_daily (normal model)
  real K600_daily_mu_mu;
  real K600_daily_mu_sigma;
  real<lower=0> K600_daily_sigma_scale;
  
  // Error distributions
  real<lower=0> err_obs_iid_sigma_scale;
  
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
  vector<lower=0>[d] K600_daily;
  
  real K600_daily_mu;
  real<lower=0> K600_daily_sigma_scaled;
  
  real<lower=0> err_obs_iid_sigma_scaled;
}

transformed parameters {
  real<lower=0> K600_daily_sigma;
  real<lower=0> err_obs_iid_sigma;
  vector[d] GPP[n];
  vector[d] ER[n];
  vector[d] KO2[n];
  vector[d] DO_mod[n];
  
  // Rescale pooling & error distribution parameters
  K600_daily_sigma = K600_daily_sigma_scale * K600_daily_sigma_scaled;
  err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled;
  
  // Model DO time series
  // * trapezoid version
  // * observation error
  // * no process error
  // * reaeration depends on DO_mod
  
  // Calculate individual process rates
  for(i in 1:n) {
    GPP[i] = GPP_daily .* frac_GPP[i];
    ER[i] = ER_daily .* frac_ER[i];
    KO2[i] = K600_daily .* KO2_conv[i];
  }
  
  // DO model
  DO_mod[1] = DO_obs_1;
  for(i in 1:(n-1)) {
    DO_mod[i+1] =
      DO_mod[i] .*
        (2.0 - KO2[i] * timestep) ./ (2.0 + KO2[i+1] * timestep) + (
        (GPP[i] + ER[i]) ./ depth[i] +
        (GPP[i+1] + ER[i+1]) ./ depth[i+1] +
        KO2[i] .* DO_sat[i] + KO2[i+1] .* DO_sat[i+1]
      ) .* (timestep ./ (2.0 + KO2[i+1] * timestep));
  }
}

model {
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

  // Hierarchical constraints on K600_daily (normal model)
  K600_daily_pred ~ normal(K600_daily_mu_mu, K600_daily_mu_sigma);
  K600_daily_sigma_scaled ~ cauchy(0, 1);
}
