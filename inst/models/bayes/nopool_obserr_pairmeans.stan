// models/bayes/stan/nopool_obserr_pairmeans.txt

data {
  
  # hyperparameters
  real GPP.daily.mu;
  real GPP.daily.sigma;
  real ER.daily.mu;
  real ER.daily.sigma;
  real K600.daily.mu;
  real K600.daily.sigma;
  
  real err.obs.sigma.min;
  real err.obs.sigma.max;
  
  int <lower=0> n; // number of observations in the day
  
  vector [n] DO.obs;
  vector [n] frac.GPP; // fractions, summing to 1, to partition GPP.daily into per-timestep rates
  vector [n] frac.ER; // fractions, summing to 1, to partition ER.daily into per-timestep rates
  vector [n] frac.D; // fractions, summing to 1, to partition K600.daily into per-timestep rates
  vector [n] depth;
  vector [n] KO2.conv;
  
}

parameters {
  
  real GPP.daily;
  real ER.daily;
  real K600.daily;
  
  // bounded prior (implied uniform) on the sd of the errors between modeled and observed DO
  real <lower=err.obs.sigma.min, upper=err.obs.sigma.max> err.obs.sigma;
  
}

model {
  
  // Declare temporary variables
  vector [n] DO.mod;
  
  // Convert daily rates to per-observation rates
  for(i in 2:n) {
    GPP[i] <- GPP.daily * (frac.GPP[i]+frac.GPP[i-1])/2 / ((depth[i]+depth[i-1])/2);
    ER[i] <- ER.daily * (frac.ER[i]+frac.ER[i-1])/2 / ((depth[i]+depth[i-1])/2);
    K[i] <- K600.daily * (KO2.conv[i]+KO2.conv[i-1])/2 * (frac.D[i]+frac.D[i-1])/2;
  }
  
  // Model DO time series
  DO.mod[1] <- DO.obs[1];
  for(i in 2:n) {
    DO.mod[i] <- (
      DO.mod[i-1] +
        GPP[i] + 
        ER[i] + 
        K[i] * (DO.sat[i] + DO.sat[i-1] - DO.mod[i-1])/2
    ) / (1 + K[i]/2);
  }
  
  // Observation error: Compare all the DO predictions to their observations
  for (i in 2:n) {
    DO.obs[i] ~ normal(DO.mod[i], err.obs.sigma);
  }
  // Prior on the sd of the errors between modeled and observed DO; this is actually unnecessary to specify
  err.obs.sigma ~ uniform(err.obs.sigma.min, err.obs.sigma.max);
  
  // Daily mean values of GPP and ER (gO2 m^-2 d^-1) and K600 (m d^-1)
  GPP.daily ~ normal(GPP.daily.mu, GPP.daily.sigma);
  ER.daily ~ normal(ER.daily.mu, ER.daily.sigma);
  K600.daily ~ normal(K600.daily.mu, K600.daily.sigma);
  
}
