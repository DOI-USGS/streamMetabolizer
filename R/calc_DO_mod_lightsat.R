#' 
#' Alison: the test does not work. Maybe I am wrong, but I think the error could bein:
#' - Line 22: in Fort Collins we wrote light=300, but shouldn't it be light=rep(300,n)??
#' - Lines 41 and 64: I tried light, light[i] and something similar to what you used with temperature (temperature=temp.water): par=light
#' 
#' Calculate a time series of DO concentrations from PMAX, ALPHA, ER, K600, and other
#' inputs
#' 
#' \code{calc_DO_mod} simulates DO with no observation or process error.
#' 
#' Accepts PMAX, ALPHA, ER, etc. and returns DO.mod. Used in many functions, including 
#' metab_mle() and predict_DO.metab_mle()
#' 
#' @inheritParams calc_DO_mod_w_fixed_error
#' @param PMAX.daily Maximum photosynthesis rate
#' @param ALPHA.daily Initial slope of the light saturation curve
#' @param light A vector of PAR of length n.
#' @export
#' @examples
#' n = 24
#' DO <- calc_DO_mod_lightsat(PMAX.daily=0.15, ALPHA.daily=0.75, ER.daily=-9, K600.daily=2.5, 
#'   DO.sat=11, depth=1, temp.water=rep(12,n), light=rep(300,n),
#'   frac.ER=rep(1/n,n), frac.D=rep(1/n,n), 
#'   DO.mod.1=8, n=n, ODE_method="Euler")
#' \dontrun{
#' plot(DO)
#' }
calc_DO_mod_lightsat <- function(
  PMAX.daily, ALPHA.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, light, frac.ER, frac.D, DO.mod.1, n, ODE_method="pairmeans", ...) {
  
 
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.mod.1
  switch(
    ODE_method,
    "Euler"={
      for(i in seq_len(n)[-1]) { # seq_len(n)[-1] is robust to n=0 and n=1, whereas 2:n is not
        DO.mod[i] <- 
          DO.mod[i-1] +
          PMAX.daily * tanh(ALPHA.daily*light[i]/PMAX.daily)/depth + 
          ER.daily * frac.ER / depth +
          convert_k600_to_kGAS(K600.daily, temperature=temp.water, gas="O2") * frac.D * (DO.sat[i] - DO.mod[i-1]) 
        }
    },
    "pairmeans"={
      # recall that DO.sat[i] is already mean(DO.sat[i], DO.sat[i-1]), and same 
      # for frac.GPP, frac.ER, frac.D, depth, and temp.water
      if(ODE_method == "pairmeans") {
        pairmean <- function(vec) {
          if(length(vec)==1) return(vec)
          (c(NA, vec[1:(length(vec)-1)]) + vec)/2
        }
        light <- pairmean(light)
        frac.ER <- pairmean(frac.ER)
        frac.D <- pairmean(frac.D)
        depth <- pairmean(depth)
        temp.water <- pairmean(temp.water)
        DO.sat <- pairmean(DO.sat)
      }
      for(i in seq_len(n)[-1]) {
        DO.mod[i] <- (
          DO.mod[i-1] +
            PMAX.daily * tanh(ALPHA.daily*light[i]/PMAX.daily)/depth + 
            ER.daily * frac.ER / depth +
            convert_k600_to_kGAS(K600.daily, temperature=temp.water, gas="O2") * frac.D * (DO.sat[i] - DO.mod[i-1]/2) 
        ) / (1 + K[i]/2)
      }
    }
  )
  
  
  # Return
  DO.mod
}
  
  
 