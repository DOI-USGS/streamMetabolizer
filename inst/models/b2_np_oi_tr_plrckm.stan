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
  real K600_lnorm_meanlog;
  real<lower=0> K600_lnorm_sdlog;
}

parameters {
  vector<lower=0>[n_days] GPP;
  vector[n_days] ER;
  vector<lower=0>[n_days] k600;
  real<lower=0> sigma;
}

transformed parameters {
  array[n_obs] vector[n_days] metab;
  array[n_obs] vector[n_days] KO2;
  for (i in 1:n_obs){
    for (t in 1:n_days){
      KO2[i,t] = (k600[t] / depth[i,t]) / ((600 / (1800.6 - (temp_water[i,t] * 120.1) + (3.7818 * temp_water[i,t]^2) - (0.047608 * temp_water[i,t]^3)))^-0.5);
    }
    metab[i] = (DO_obs_up[i] + GPP .* light[i] ./ depth[i] + ER .* travel_time[i] ./ depth[i] + (KO2[i] .* travel_time[i] .* (DO_sat_up[i] - DO_obs_up[i] + DO_sat_down[i]) / 2)) ./ (1 + (KO2[i] .* travel_time[i]) / 2);
  }
}

model {
  for (i in 1:n_obs){
    DO_obs_down[i] ~ normal(metab[i], sigma);
  }
  for (i in 1:n_days){
    k600[i] ~ lognormal(K600_lnorm_meanlog, K600_lnorm_sdlog);
  }
}
