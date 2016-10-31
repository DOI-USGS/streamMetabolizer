// b_np_pc_tr_plrcko.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  // Error distributions
  real err_proc_acor_phi_alpha;
  real err_proc_acor_phi_beta;
  real<lower=0> err_proc_acor_sigma_scale;
  
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
  
  real<lower=0, upper=1> err_proc_acor_phi;
  real<lower=0> err_proc_acor_sigma_scaled;
  
  vector[d] err_proc_acor_inc[n+1];
}

transformed parameters {
  vector[d] DO_mod_partial_sigma[n];
  real<lower=0> err_proc_acor_sigma;
  vector[d] GPP[n];
  vector[d] ER[n];
  vector[d] KO2[n];
  vector[d] DO_mod[n];
  vector[d] err_proc_acor[n];
  
  // Rescale pooling & error distribution parameters
  err_proc_acor_sigma = err_proc_acor_sigma_scale * err_proc_acor_sigma_scaled;
  
  // Model DO time series
  // * trapezoid version
  // * no observation error
  // * autocorrelated process error
  // * reaeration depends on DO_obs
  
  err_proc_acor[1] = err_proc_acor_inc[1];
  for(i in 1:(n-1)) {
    err_proc_acor[i+1] = err_proc_acor_phi * err_proc_acor[i] + err_proc_acor_inc[i+1];
  }
  
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
      DO_obs[i] + (
        - KO2[i] .* DO_obs[i] - KO2[i+1] .* DO_obs[i+1] +
        (GPP[i] + ER[i] + err_proc_acor[i]) ./ depth[i] +
        (GPP[i+1] + ER[i+1] + err_proc_acor[i+1]) ./ depth[i+1] +
        KO2[i] .* DO_sat[i] + KO2[i+1] .* DO_sat[i+1]
      ) * (timestep / 2.0);
  }
}

model {
  // Process error
  for(i in 2:n) {
    // Autocorrelated process error
    err_proc_acor_inc[i-1] ~ normal(0, err_proc_acor_sigma);
  }
  // Autocorrelation (phi) & SD (sigma) of the process errors
  err_proc_acor_phi ~ beta(err_proc_acor_phi_alpha, err_proc_acor_phi_beta);
  err_proc_acor_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
}
