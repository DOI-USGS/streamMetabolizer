#' \code{calc_DO_mod_w_fixed_error} simulates DO with vectors of observation 
#' and/or process error.
#' 
#' @rdname calc_DO_mod
#'   
#' @param GPP.daily One GPP rate per day (mg O2 / L / d)
#' @param ER.daily One ER rate per day (mg O2 / L / d), always non-positive.
#' @param K600.daily One K600 per day (1 / d)
#' @param DO.sat dissolved oxygen concentrations if the water were at 
#'   equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}. Calculate using 
#'   \link{calc_DO_at_sat}
#' @param depth stream depth, \eqn{m}{m}.
#' @param temp.water stream temperature in degC
#' @param frac.GPP the fraction of daily GPP to apply to each timestep
#' @param frac.ER the fraction of daily ER to apply to each timestep
#' @param frac.D the fraction of daily D to apply to each timestep
#' @param DO.mod.1 the first DO.obs value, to which the first DO.mod value will 
#'   be set
#' @param n number of DO.mod values to produce
#' @param ... additional arguments passed to other variants on calc_DO_mod
#' @param err.obs A vector of observation errors of length n. Observation errors
#'   are those applied to DO.mod after generating the full time series of 
#'   modeled values.
#' @param err.proc A vector process errors of length n. Process errors are 
#'   applied at each time step, and therefore propagate into the next timestep.
#' @param ODE_method character specifying the numerical integration method to 
#'   use. The default is pairmeans, where the change in DO between times t-1 and
#'   t is a function of the mean DO.sat between times t-1 and t, mean temp.water
#'   between times t-1 and t, mean light between times t-1 and t, etc. The Euler
#'   method is currently more common in the literature, with each time step 
#'   depending entirely on DO.sat, GPP, etc. at time t. Both methods are
#'   imprecise and fast, relative to Runge-Kutta or other numerical integration
#'   methods
#' @export
#' @examples
#' fr <- rep(1/48,48) # shorthand for readability of next lines
#' gfr <- sin((1:48)/(4*pi))^8; gfr <- gfr/sum(gfr)
#' obsrnorm <- stats::rnorm(48,0,0.4)
#' procrnorm <- stats::rnorm(48,0,0.2)
#' preds <- lapply(c("Euler","pairmeans"), function(method) {
#'  data.frame(
#'   t=1:48,
#'   method=method,
#'   no_err = calc_DO_mod_w_fixed_error(
#'    10, -13, 2.5, 14, 1, rep(12,48), gfr, fr, fr, 11, 48, 0, 0, method),
#'   obs_err_const = calc_DO_mod_w_fixed_error(
#'    10, -13, 2.5, 14, 1, rep(12,48), gfr, fr, fr, 11, 48, 0.1, 0, method),
#'   obs_err_rnorm = calc_DO_mod_w_fixed_error(
#'    10, -13, 2.5, 14, 1, rep(12,48), gfr, fr, fr, 11, 48, obsrnorm, 0, method),
#'   proc_err_const = calc_DO_mod_w_fixed_error(
#'    10, -13, 2.5, 14, 1, rep(12,48), gfr, fr, fr, 11, 48, 0, 0.01, method),
#'   proc_err_rnorm = calc_DO_mod_w_fixed_error(
#'    10, -13, 2.5, 14, 1, rep(12,48), gfr, fr, fr, 11, 48, 0, procrnorm, method),
#'   all_err = calc_DO_mod_w_fixed_error(
#'    10, -13, 2.5, 14, 1, rep(12,48), gfr, fr, fr, 11, 48, 
#'    obsrnorm, procrnorm, method),
#'   stringsAsFactors=FALSE
#'  )
#' })
#' \dontrun{
#' library(ggplot2); library(tidyr); library(dplyr)
#' ggplot(gather(bind_rows(preds), type, value, no_err, obs_err_const, obs_err_rnorm, 
#'               proc_err_const, proc_err_rnorm, all_err), 
#'   aes(x=t, y=value, color=type)) + geom_line(aes(linetype=method)) + theme_bw() +
#'   facet_wrap(~ type)
#' ggplot(preds %>% bind_rows %>% 
#'   mutate(perfect = rep(no_err[method=="pairmeans"], times=2)) %>%
#'   gather(type, value, no_err, obs_err_const, obs_err_rnorm, 
#'          proc_err_const, proc_err_rnorm, all_err), 
#'   aes(x=t, y=value-perfect, color=type)) + geom_line(aes(linetype=method)) + 
#'   theme_bw() + facet_wrap(~ type)
#' }
calc_DO_mod_w_fixed_error <- function(
  GPP.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n,
  err.obs=rep(0, n), err.proc=rep(0, n), ODE_method=c("pairmeans","Euler"), ...) {
  
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
  
  # make sure anything required in the core modeling loop has n observations
  if(length(DO.sat) != n) DO.sat <- rep(DO.sat, length.out=n) # this would be odd
  if(length(err.proc) != n) err.proc <- rep(err.proc, length.out=n) # this is more probable. err.obs can stay length 1 if it is
  
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.mod.1
  switch(
    ODE_method,
    "Euler"={
      for(i in seq_len(n)[-1]) { # seq_len(n)[-1] is robust to n=0 and n=1, whereas 2:n is not
        DO.mod[i] <- 
          DO.mod[i-1] +
          GPP[i] + 
          ER[i] + 
          K[i] * (DO.sat[i] - DO.mod[i-1]) +
          err.proc[i]
      }
    },
    "pairmeans"={
      # recall that DO.sat[i] is already mean(DO.sat[i], DO.sat[i-1]), and same 
      # for frac.GPP, frac.ER, frac.D, depth, and temp.water
      for(i in seq_len(n)[-1]) {
        DO.mod[i] <- (
          DO.mod[i-1] +
          GPP[i] + 
          ER[i] + 
          K[i] * (DO.sat[i] - DO.mod[i-1]/2) +
          err.proc[i]
        ) / (1 + K[i]/2)
      }
    }
  )
  DO.mod <- DO.mod + err.obs
  
  # Return
  DO.mod
}