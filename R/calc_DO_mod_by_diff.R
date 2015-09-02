#' Simulates DO using a process error model, i.e., AKA one-step-ahead prediction
#' 
#' Accepts GPP, ER, etc. and returns DO.mod. This is a process model, in the 
#' sense that DO.mod[t] is a function of DO.obs[t-1], so each DO.mod value 
#' essentially describes only the processes occurring between times t-1 and t, 
#' and the difference between DO.mod[t] and DO.obs[t] is interpreted as error in
#' the process (as a rate of mg/L per timestep)
#' 
#' @inheritParams calc_DO_mod
#' @param DO.obs a vector of observed DO, from which values at time t will be 
#'   used to model values at time t+1
#' @export
calc_DO_mod_by_diff <- function(
  GPP.daily, ER.daily, K600.daily, DO.obs, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, n, ...) {
  
  # partition GPP and ER into their timestep-specific rates (mg/L/timestep at
  # each timestep)
  GPP <- GPP.daily * frac.GPP / depth
  ER <- ER.daily * frac.ER / depth
  K <- convert_k600_to_kGAS(K600.daily, temperature=temp.water, gas="O2") * frac.D
  
  # make sure anything in the following loop has n observations
  if(length(DO.sat) != n) DO.sat <- rep(DO.sat, length.out=n) # this would be odd
  
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.obs[1]
  for(i in 2:n) {
    DO.mod[i] <- DO.obs[i-1] +
      GPP[i] + 
      ER[i] + 
      K[i] * (DO.sat[i] - DO.obs[i-1])
  } 
  
  # Return
  DO.mod
}