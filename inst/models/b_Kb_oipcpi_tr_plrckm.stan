// b_Kb_oipcpi_tr_plrckm.stan

data {
  // Parameters of priors on metabolism
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  
  // Parameters of hierarchical priors on K600_daily (binned model)
  int <lower=1> b; # number of K600_lnQ_nodes
  real K600_lnQ_nodediffs_sdlog;
  vector[b] K600_lnQ_nodes_meanlog;
  vector[b] K600_lnQ_nodes_sdlog;
  real<lower=0> K600_daily_sdlog_scale;
  
  // Error distributions
  real<lower=0> err_obs_iid_sigma_scale;
  real err_proc_acor_phi_alpha;
  real err_proc_acor_phi_beta;
  real<lower=0> err_proc_acor_sigma_scale;
  real<lower=0> err_proc_iid_sigma_scale;
  
  // Data dimensions
  int<lower=1> d; # number of dates
  int<lower=1> n; # number of observations per date
  
  // Daily data
  vector[d] DO_obs_1;
  int<lower=1,upper=b> lnQ_bins[2,d];
  vector<lower=0,upper=1>[d] lnQ_bin_weights[2];
  
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
  
  vector[b] lnK600_lnQ_nodes;
  real<lower=0> K600_daily_sdlog_scaled;
  
  real<lower=0> err_obs_iid_sigma_scaled;
  real<lower=0, upper=1> err_proc_acor_phi;
  real<lower=0> err_proc_acor_sigma_scaled;
  real<lower=0> err_proc_iid_sigma_scaled;
  
  vector[d] err_proc_iid[n-1];
  vector[d] err_proc_acor_inc[n+1];
  vector[d] DO_mod[n];
}

transformed parameters {
  real<lower=0> K600_daily_sdlog;
  vector[d] K600_daily_pred;
  real<lower=0> err_obs_iid_sigma;
  vector[d] DO_mod_partial_sigma[n];
  real<lower=0> err_proc_acor_sigma;
  real<lower=0> err_proc_iid_sigma;
  vector[d] GPP[n];
  vector[d] ER[n];
  vector[d] KO2[n];
  vector[d] DO_mod_partial[n];
  vector[d] err_proc_acor[n];
  
  // Rescale pooling & error distribution parameters
  K600_daily_sdlog = K600_daily_sdlog_scale * K600_daily_sdlog_scaled;
  err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled;
  err_proc_acor_sigma = err_proc_acor_sigma_scale * err_proc_acor_sigma_scaled;
  err_proc_iid_sigma = err_proc_iid_sigma_scale * err_proc_iid_sigma_scaled;
  
  // Hierarchical, binned model of K600_daily
  K600_daily_pred = exp(lnK600_lnQ_nodes[lnQ_bins[1]] .* lnQ_bin_weights[1] + 
                        lnK600_lnQ_nodes[lnQ_bins[2]] .* lnQ_bin_weights[2]);
  
  // Model DO time series
  // * trapezoid version
  // * observation error
  // * IID and autocorrelated process error
  // * reaeration depends on DO_mod
  
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
  for(i in 1:(n-1)) {
    DO_mod_partial[i+1] =
      DO_mod[i] .*
        (2.0 - KO2[i] * timestep) ./ (2.0 + KO2[i+1] * timestep) + (
        (GPP[i] + ER[i] + err_proc_acor[i]) ./ depth[i] +
        (GPP[i+1] + ER[i+1] + err_proc_acor[i+1]) ./ depth[i+1] +
        KO2[i] .* DO_sat[i] + KO2[i+1] .* DO_sat[i+1]
      ) .* (timestep ./ (2.0 + KO2[i+1] * timestep));
    for(j in 1:d) {
      DO_mod_partial_sigma[i+1,j] = err_proc_iid_sigma * 
        sqrt(pow(depth[i,j], -2) + pow(depth[i+1,j], -2)) .*
        (timestep / (2.0 + KO2[i+1,j] * timestep));
    }
  }
}

model {
  // Process error
  for(i in 2:n) {
    // Independent, identically distributed process error
    DO_mod[i] ~ normal(DO_mod_partial[i], DO_mod_partial_sigma[i]);
    // Autocorrelated process error
    err_proc_acor_inc[i-1] ~ normal(0, err_proc_acor_sigma);
  }
  // SD (sigma) of the IID process errors
  err_proc_iid_sigma_scaled ~ cauchy(0, 1);
  // Autocorrelation (phi) & SD (sigma) of the process errors
  err_proc_acor_phi ~ beta(err_proc_acor_phi_alpha, err_proc_acor_phi_beta);
  err_proc_acor_sigma_scaled ~ cauchy(0, 1);
  
  // Independent, identically distributed observation error
  for(i in 2:n) {
    DO_obs[i] ~ normal(DO_mod[i], err_obs_iid_sigma);
  }
  // SD (sigma) of the observation errors
  err_obs_iid_sigma_scaled ~ cauchy(0, 1);
  
  // Daily metabolism priors
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ lognormal(K600_daily_pred, K600_daily_sdlog);

  // Hierarchical constraints on K600_daily (binned model)
  lnK600_lnQ_nodes ~ normal(K600_lnQ_nodes_meanlog, K600_lnQ_nodes_sdlog);
  for(k in 2:b) {
    lnK600_lnQ_nodes[k] ~ normal(lnK600_lnQ_nodes[k-1], K600_lnQ_nodediffs_sdlog);
  }
  K600_daily_sdlog_scaled ~ cauchy(0, 1);
}
