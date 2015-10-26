// This version follows from nopool_pcpi_Euler_b2.stan but attempts to 
// gain speed by vectorization, transforming data & parameters early, etc.

data {
 
  real GPP_daily_mu;
  real GPP_daily_sigma;
  real ER_daily_mu;
  real ER_daily_sigma;
  real K600_daily_mu;
  real K600_daily_sigma;
  
  real err_proc_acor_phi_min;
  real err_proc_acor_phi_max;
  real err_proc_acor_sigma_min;
  real err_proc_acor_sigma_max;
  real err_proc_iid_sigma_min;
  real err_proc_iid_sigma_max;
  
  int <lower=0> n;

  vector [n] DO_obs;
  vector [n] DO_sat;
  vector [n] frac_GPP;
  vector [n] frac_ER;
  vector [n] frac_D;
  vector [n] depth;
  vector [n] KO2_conv;
  
}

// manual says: the statements in the transformed data block are only ever evaluated once
transformed data {
  
  vector [n-1] coef_GPP;
  vector [n-1] coef_ER;
  vector [n-1] coef_K600;
  vector [n-1] dDOdt_obs;
  
  // coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)
  for(i in 1:(n-1)) {
    coef_GPP[i]  <- frac_GPP[i] / depth[i];
    coef_ER[i]   <- frac_ER[ i] / depth[i];
    coef_K600[i] <- KO2_conv[i] * frac_D[i] * 
      (DO_sat[i] - DO_obs[i]);
    dDOdt_obs[i] <- DO_obs[i+1] - DO_obs[i];
  }  
  
}

parameters {

  real GPP_daily;
  real ER_daily;
  real K600_daily;
  real <lower=err_proc_acor_phi_min,   upper=err_proc_acor_phi_max>   err_proc_acor_phi;
  real <lower=err_proc_acor_sigma_min, upper=err_proc_acor_sigma_max> err_proc_acor_sigma;
  real <lower=err_proc_iid_sigma_min,  upper=err_proc_iid_sigma_max>  err_proc_iid_sigma;
  vector[n-1] err_proc_acor;

}

// manual says: evaluated once per leapfrog step
transformed parameters {
  
  vector [n-1] dDOdt_mod;
  
  // DO process model. Can be vectorized & pre-transformed b/c we're 
  // using DO_obs rather than DO_mod to compute the DO deficit
  dDOdt_mod <- 
    err_proc_acor + 
    GPP_daily * coef_GPP + 
    ER_daily * coef_ER +
    K600_daily * coef_K600;
    
}

model {
  
  // Autocorrelated process error (something seems to be implied for err_proc_acor[1])
  // err_proc_acor[1] ~ normal(0, err_proc_acor_sigma);
  for(i in 2:(n-1)) {
    err_proc_acor[i] ~ normal(err_proc_acor_phi*err_proc_acor[i-1], err_proc_acor_sigma);
  }
  // Prior on the autocorrelation & sd of the process errors (uniforms implied in parameters declaration)
  //err_proc_acor_phi ~ uniform(err_proc_acor_phi_min, err_proc_acor_phi_max);
  //err_proc_acor_sigma ~ uniform(err_proc_acor_sigma_min, err_proc_acor_sigma_max);
  
  // Uncorrelated (IID) process error: Compare all the dDOdt predictions to their observations
  dDOdt_obs ~ normal(dDOdt_mod, err_proc_iid_sigma);
  // Prior on the sd of the errors between modeled and observed DO (uniform implied in parameters declaration)
  //err_proc_iid_sigma ~ uniform(err_proc_iid_sigma_min, err_proc_iid_sigma_max);
  
  // Daily mean values of GPP and ER (gO2 m^-2 d^-1) and K600 (m d^-1)
  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma);
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma);
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma);
  
}