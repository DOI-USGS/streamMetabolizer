#' \code{calc_DO_mod_by_diff} simulates DO using a process error model, i.e., 
#' AKA one-step-ahead prediction. In this model, DO.mod[t] is a function of 
#' DO.obs[t-1], so each DO.mod value essentially describes only the processes 
#' occurring between times t-1 and t. The difference between DO.mod[t] and 
#' DO.obs[t] can therefore be interpreted as error in the process estimate (as a
#' rate of mg/L per timestep)
#' 
#' @rdname calc_DO_mod
#'   
#' @inheritParams calc_DO_mod_w_fixed_error
#' @param DO.obs a vector of observed DO, from which values at time t will be 
#'   used to model values at time t+1
#' @export
#' @examples 
#' fr <- rep(1/48,48) # shorthand for readability of next lines
#' gfr <- sin((1:48)/(4*pi))^8; gfr <- gfr/sum(gfr)
#' do <- 2*sin((1:48 - 10)/(4*pi))^4 + 10
#' domod_pm <- calc_DO_mod_by_diff(10, -13, 2.5, do, 14, 1, rep(12,48), gfr, fr, fr, 48, "pairmeans")
#' domod_eu <- calc_DO_mod_by_diff(10, -13, 2.5, do, 14, 1, rep(12,48), gfr, fr, fr, 48, "Euler")
#' \dontrun{
#' plot(domod_eu)
#' points(domod_pm, col='red')
#' }
calc_DO_mod_by_diff <- function(
  GPP.daily, ER.daily, K600.daily, DO.obs, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, n, ODE_method=c("pairmeans","Euler"), ...) {
  
  # if requested, use the means of i and i-1 for DO.sat, DO.mod, temperature, and light
  ODE_method <- match.arg(ODE_method)
  if(ODE_method == "pairmeans") {
    pairmean <- function(vec) {
      if(length(vec)==1) return(vec)
      (c(NA, vec[1:(length(vec)-1)]) + vec)/2
    }
    frac.GPP <- pairmean(frac.GPP)
    frac.ER <- pairmean(frac.ER)
    frac.D <- pairmean(frac.D)
    depth <- pairmean(depth)
    temp.water <- pairmean(temp.water)
    DO.sat <- pairmean(DO.sat)
  }
  
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
  switch(
    ODE_method,
    "Euler"={
      for(i in seq_len(n)[-1]) {
        DO.mod[i] <- DO.obs[i-1] +
          GPP[i] + 
          ER[i] + 
          K[i] * (DO.sat[i] - DO.obs[i-1])
      }
    },
    "pairmeans"={
      # recall that DO.sat[i] is already mean(DO.sat[i], DO.sat[i-1]), and same 
      # for frac.GPP, frac.ER, frac.D, depth, and temp.water
      for(i in seq_len(n)[-1]) {
        DO.mod[i] <- DO.obs[i-1] +
          GPP[i] + 
          ER[i] + 
          K[i] * (DO.sat[i] - (DO.obs[i] + DO.obs[i-1])/2)
      }
    }
  )
  
  # Return
  DO.mod
}