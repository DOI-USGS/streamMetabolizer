// This script intended to stay as close to nopool_pcpi_Euler_b1.stan
// as possible except for renaming to streamMetabolizer standards and
// transforming data where the b1 version implies it was done before
// calling the model. Spacing in this file should correspond closely
// to that in nopool_pcpi_Euler_b1.stan.

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
  
  int <lower=0> n; // was T

  vector [n] DO_obs; // was y
  
  vector [n] DO_sat; // was oxysat
  vector [n] frac_GPP; // was light
  vector [n] frac_ER; // was 0.006944
  vector [n] frac_D; // was 0.006944
  vector [n] depth; // was z
  vector [n] KO2_conv; // was Kc
  
}

transformed data {
  
  vector [n] dDOdt_obs; // was dodt; was calculated before model call
  
  for(i in 2:n) {
    dDOdt_obs[i-1] <- (DO_obs[i] - DO_obs[i-1]);
  }  
  dDOdt_obs[n] <- dDOdt_obs[n-1];
  
}

parameters {
  
  real GPP_daily; // was GPP
  real ER_daily; // was ER
  real K600_daily; // was K
  real <lower=err_proc_acor_phi_min, upper=err_proc_acor_phi_max> err_proc_acor_phi; // was phi<0,1>
  real <lower=err_proc_iid_sigma_min, upper=err_proc_iid_sigma_max> err_proc_iid_sigma; // was sigobs<0,>
  real <lower=err_proc_acor_sigma_min, upper=err_proc_acor_sigma_max> err_proc_acor_sigma; // was sigproc<0,>
  vector[n] err_proc_acor; // was eta
  
}

model {

  vector [n] dDOdt_mod; // was mdodt
  
  for (i in 1:n){
    dDOdt_mod[i] <- 
      err_proc_acor[i] + 
      ((GPP_daily/depth[i])*(frac_GPP[i]/sum(frac_GPP))) + 
      (ER_daily*frac_ER[i]/depth[i]) +
      (KO2_conv[i]*K600_daily * frac_D[i]*(DO_sat[i]-DO_obs[i])); 
  }

  GPP_daily ~ normal(GPP_daily_mu, GPP_daily_sigma); // was normal(0, 10)
  ER_daily ~ normal(ER_daily_mu, ER_daily_sigma); // was normal(-10, 10)
  K600_daily ~ normal(K600_daily_mu, K600_daily_sigma); // was normal(10, 10)
  
  for (i in 1:n){
    dDOdt_obs[i] ~ normal(dDOdt_mod[i], err_proc_iid_sigma);
  }           

  for (i in 2:n){
    err_proc_acor[i] ~ normal(err_proc_acor_phi*err_proc_acor[i-1], err_proc_acor_sigma);
  }
  
}