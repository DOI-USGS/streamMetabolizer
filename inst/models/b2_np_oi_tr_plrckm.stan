// b2_np_oi_tr_plrckm.stan

data {
  int<lower=1> n_obs;  // number of total do observations
  int<lower=1> n_days; // number of days
  array[n_obs] vector[n_days] DO_obs_up;
  array[n_obs] vector[n_days] DO_sat_up;
  array[n_obs] vector[n_days] DO_obs_down;
  array[n_obs] vector[n_days] DO_sat_down;
  array[n_obs] vector[n_days] light;
  array[n_obs] vector[n_days] depth;
  array[n_obs] vector[n_days] temp_water;
  array[n_obs] vector[n_days] travel_time;
  real GPP_daily_mu;
  real<lower=0> GPP_daily_sigma;
  real ER_daily_mu;
  real<lower=0> ER_daily_sigma;
  real k600_velocity_lnorm_meanlog;
  real<lower=0> k600_velocity_lnorm_sdlog;
}

parameters {
  vector<lower=0>[n_days] GPP_daily;
  vector[n_days] ER_daily;
  vector<lower=0>[n_days] k600_velocity; // gas exchange velocity, m d^-1, constant within a day
  real<lower=0> sigma;
}

transformed parameters {
  array[n_obs] vector[n_days] K600; // rate-scale gas exchange, d^-1
  array[n_obs] vector[n_days] metab;
  array[n_obs] vector[n_days] KO2;
  for (i in 1:n_obs){
    for (t in 1:n_days){
      K600[i,t] = k600_velocity[t] / depth[i,t];
      KO2[i,t] = K600[i,t] / ((600 / (1800.6 - (temp_water[i,t] * 120.1) + (3.7818 * temp_water[i,t]^2) - (0.047608 * temp_water[i,t]^3)))^-0.5);
    }
    metab[i] = (DO_obs_up[i] + GPP_daily .* light[i] ./ depth[i] + ER_daily .* travel_time[i] ./ depth[i] + (KO2[i] .* travel_time[i] .* (DO_sat_up[i] - DO_obs_up[i] + DO_sat_down[i]) / 2)) ./ (1 + (KO2[i] .* travel_time[i]) / 2);
  }
}

model {
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  for (i in 1:n_days){
    k600_velocity[i] ~ lognormal(k600_velocity_lnorm_meanlog, k600_velocity_lnorm_sdlog);
  }
  for (i in 1:n_obs){
    DO_obs_down[i] ~ normal(metab[i], sigma);
  }
}

generated quantities {
  // daily mean rate-scale K600, d^-1: divide the day's velocity-scale
  // k600 by each timestep's depth first, then average over the day
  vector[n_days] K600_daily;
  for (t in 1:n_days){
    real K600_sum = 0;
    for (i in 1:n_obs){
      K600_sum += K600[i,t];
    }
    K600_daily[t] = K600_sum / n_obs;
  }
}
