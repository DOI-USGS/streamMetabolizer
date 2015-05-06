#' @title Basic metabolism model fitting function
#' @description Fits a model to estimate GPP and ER from input data on DO, 
#'   temperature, light, etc.
#'   
#' @param data data.frame with columns \itemize{
#'   
#'   \item{ \code{date.time} date-time values in POSIXct format
#'   
#'   \item{ \code{DO.obs} dissolved oxygen concentration observations, \eqn{mg 
#'   O[2] L^{-1}}{mg O2 / L}}
#'   
#'   \item{ \code{DO.sat} dissolved oxygen concentrations if the water were at 
#'   equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}}. Calculate using 
#'   \link{calc_DO_at_sat}}
#'   
#'   \item{ \code{depth} stream depth, \eqn{m}{m}}.
#'   
#'   \item{ \code{k.O2} gas exchange coefficients, \eqn{m d^{-1}}{m / d}}.
#'   
#'   \item{ \code{light} photosynthetically active radiation, \eqn{\mu mol\ 
#'   m^{-2} s^{-1}}{micro mols / m^2 / s}}
#'   
#'   }
#' @param ... additional arguments
#' @return A metab_simple object containing the fitted model.
#'   
#' @author Alison Appling, Jordan Read; modeled on lakeMetabolizer
#' @examples
#' \dontrun{
#'  metab_simple(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_simple <- function(data, ...) {
  
  # Check data for correct column names
  expected.colnames <- c("date.time","DO.obs","DO.sat","depth","k.O2","light")
  if(!all(expected.colnames %in% names(data))) {
    stop(paste0("data must contain (at least) columns with the names ", paste0(expected.colnames, collapse=", ")))
  }
  
  # Require that date.time have a single, consistent time step which we'll extract
  timestep.days <- mean(as.numeric(diff(data$date.time), units="days"))
  timestep.deviations <- diff(range(as.numeric(diff(data$date.time), units="days")))
  if((timestep.deviations / timestep.days) > 0.001) { # threshold: max-min timestep length can't be more than 1% of mean timestep length
    stop("expected essentially one timestep length within date.time")
  }
  if(abs(timestep.days * nrow(data) - 1) > 0.001) { # threshold: all timesteps per day must add up to 1, plus or minus 0.001
    stop("number * duration of time intervals doesn't sum to 1 day")
  }
  
  #' Return the likelihood value for a given set of parameters and observations
  #' 
  #' From ?nlm, this function should be "the function to be minimized, returning
  #' a single numeric value. This should be a function with first argument a 
  #' vector of the length of p followed by any other arguments specified by the 
  #' ... argument."
  #' 
  #' @param params a vector of length 2, where the first element is GPP and the second element is ER (both mg/L/d)
  onestation_negloglik <- function(params, DO.obs, DO.sat, depth, k.O2, frac.GPP, frac.ER, frac.D) {
    
    # Count how many DO observations/predictions there should be
    n <- length(DO.obs)
    
    # Parse params vector (passed from nlm) and produce DO.mod estimates
    DO.mod <- .core_model_metab_simple(
      GPP.daily=params[1], ER.daily=params[2], 
      DO.sat, depth, k.O2, frac.GPP, frac.ER, frac.D, DO.obs[1], n)
    
    # calculate & return the negative log likelihood of DO.mod values relative
    # to DO.obs values. equivalent to Bob's original code & formula at
    # http://www.statlect.com/normal_distribution_maximum_likelihood.htm
    diffs.sq <- (DO.obs-DO.mod)^2 
    sigma.sq <- sum(diffs.sq)/n
    (n/2)*log(sigma.sq) + (n/2)*log(2*pi) + (1/(2*sigma.sq))*sum(diffs.sq)
  }
  
  # Calculate metabolism by non linear minimization of an MLE function
  river.mle <- with(
    data, 
    nlm(onestation_negloglik, p=c(GPP=3,ER=-5), 
        DO.obs=DO.obs, DO.sat=DO.sat, depth=depth, k.O2=k.O2, 
        frac.GPP=light/sum(light), frac.ER=timestep.days, frac.D=timestep.days)
  )
  
  # Package and return results. There are no args, just data, for this
  # particular model
  new("metab_simple", 
      fit=river.mle,
      args=list(),
      data=data[expected.colnames],
      pkg_version=as.character(packageVersion("streamMetabolizer")))
}

#' This is the workhorse for metab_simple, used in both metab_simple() and 
#' predict_DO.metab_simple(). It accepts GPP, ER, etc. and returns DO.mod.
#' 
#' @param GPP.daily One GPP rate per day (mg O2 / L / d)
#' @param ER.daily One ER rate per day (mg O2 / L / d), always non-positive.
#' @param date.time date-time values in POSIXct format
#' @param DO.obs dissolved oxygen concentration observations, \eqn{mg O[2]
#'   L^{-1}}{mg O2 / L}
#' @param DO.sat dissolved oxygen concentrations if the water were at 
#'   equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}. Calculate using 
#'   \link{calc_DO_at_sat}
#' @param depth stream depth, \eqn{m}{m}.
#' @param k.O2 gas exchange coefficients, \eqn{m d^{-1}}{m / d}.
#' @param light photosynthetically active radiation, \eqn{\mu mol\ m^{-2}
#'   s^{-1}}{micro mols / m^2 / s}
#' @param DO.obs.1 the first DO.obs value, to which the first DO.mod value will be set
#' @param n number of DO.mod values to produce
#'   
#' @name core_model_metab_simple
#' @keywords internal
.core_model_metab_simple <- function(GPP.daily, ER.daily, DO.sat, depth, k.O2, frac.GPP, frac.ER, frac.D, DO.obs.1, n) {
  # partition GPP and ER into their timestep-specific rates (mg/L/timestep at
  # each timestep)
  GPP <- GPP.daily * frac.GPP / depth
  ER <- ER.daily * frac.ER / depth
  
  # make sure anything in the following loop has n observations
  if(length(k.O2) != n) k.O2 <- rep(k.O2, length.out=n) # this is possible but unlikely
  if(length(DO.sat) != n) DO.sat <- rep(DO.sat, length.out=n) # this would be odd
  if(length(frac.D) != n) frac.D <- rep(frac.D, length.out=n) # this is probable. guaranteed? not if there's a missing timestep.
  
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.obs.1
  for(i in 2:n) {
    DO.mod[i] <- DO.mod[i-1] +
      GPP[i] + 
      ER[i] + 
      k.O2[i] * (DO.sat[i] - DO.mod[i-1]) * frac.D[i]
  } 
  
  # Return
  DO.mod
}


#### metab_simple class ####

#' A metabolism model class specific to a very simple model.
#' 
#' metab_simple models expect K as an input and use non-linear minimization to
#' fit values of GPP and ER for a given DO curve.
#' 
#' @exportClass metab_simple
#' @family metab.model.classes
setClass(
  "metab_simple", 
  contains="metab_model"
)


#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily predictions of GPP, ER, and NEP.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @export
#' @family predict_metab
predict_metab.metab_simple <- function(metab_model) {
  
  # at the moment this will only work if there's only one date
  date <- unique(as.Date(format(get_data(metab_model)$date.time, "%Y-%m-%d")))
  GPP <- get_fit(metab_model)$estimate[1]
  ER <- get_fit(metab_model)$estimate[2]
  
  data.frame(date=date, GPP=GPP, ER=ER, NEP=GPP+ER)
}



#' Make dissolved oxygen predictions from a fitted metab_model.
#' 
#' Makes fine-scale predictions of dissolved oxygen using fitted coefficients, 
#' etc. from the metabolism model.
#' 
#' @import dplyr
#' @inheritParams predict_DO
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @export
#' @family predict_DO
predict_DO.metab_simple <- function(metab_model) {

  # get the metabolism (GPP, ER) estimates
  metab_ests <- predict_metab(metab_model)
  
  # re-process the input data with the metabolism estimates to predict dissolved
  # oxygen
  . <- date.time <- ".dplyr.var"
  get_data(metab_model) %>%
    mutate(date=as.Date(format(date.time, "%Y-%m-%d"))) %>%
    group_by(date) %>%
    do(with(., {
      # prepare auxiliary data
      n <- length(date.time)
      timestep.days <- mean(as.numeric(diff(date.time), units="days"))
      frac.GPP=light/sum(light)
      frac.ER=timestep.days
      frac.D=timestep.days
      
      # produce DO.mod estimates for today's GPP and ER
      DO.mod <- .core_model_metab_simple(
        GPP.daily=metab_ests[metab_ests$date==date[1], "GPP"], ER.daily=metab_ests[metab_ests$date==date[1], "ER"], 
        DO.sat, depth, k.O2, frac.GPP, frac.ER, frac.D, DO.obs[1], n)
      
      data.frame(., DO.mod=DO.mod)
    }))
    
}
