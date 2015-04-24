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
#'   \item{ \code{DO.deficit} dissolved oxygen saturation deficit values,
#'   positive when observed DO is less than the equilibrium DO concentration,
#'   \eqn{mg O[2] L^{-1}}{mg O2 / L}}. Calculate using \link{calc_DO_deficit}}
#'   
#'   \item{ \code{depth} stream depth, \eqn{m}{m}}.
#'   
#'   \item{ \code{k.O2} gas exchange coefficients, \eqn{m d^{-1}}{m / d}}.
#'   
#'   \item{ \code{light} photosynthetically active radiation, \eqn{\mu mol\ m^{-2}
#'   s^{-1}}{micro mols / m^2 / s}}
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
  expected.colnames <- c("date.time","DO.obs","DO.deficit","depth","k.O2","light")
  if(!all(expected.colnames %in% names(data))) {
    stop(paste0("data must contain (at least) columns with the names ", paste0(expected.colnames, collapse=", ")))
  }
  
  # Require that date.time have a single, consistent time step which we'll extract
  timestep.days <- unique(as.numeric(diff(data$date.time), units="days"))
  if(length(timestep.days) != 1) {
    stop("expected one and only one timestep length within date.time")
  }
  if(timestep.days * nrow(data) != 1) {
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
  #' @param 
  onestation_negloglik <- function(params, DO.obs, DO.deficit, depth, k.O2, frac.GPP, frac.ER, frac.D) {
    
    # Parse params vector (passed from nlm)
    GPP <- params[1]
    ER <- params[2]
    
    # Model DO with given params
    DO.mod <- numeric(length(DO.obs))
    DO.mod[1] <- DO.obs[1]
    DO.delta <- 
      GPP * frac.GPP / depth + 
      ER * frac.ER / depth + 
      k.O2 * DO.deficit * frac.D # Bob - pros & cons of using DO.mod[t-1] instead of DO.obs[t] or DO.obs[t-1] (i.e., from DO.deficit) here?
    DO.mod[-1] <- DO.mod[1] + cumsum(DO.delta[-1])
    
    # calculate & return the negative log likelihood of DO.mod values relative
    # to DO.obs values. equivalent to Bob's original code & formula at
    # http://www.statlect.com/normal_distribution_maximum_likelihood.htm
    n <- length(DO.obs)
    diffs.sq <- (DO.obs-DO.mod)^2 
    sigma.sq <- sum(diffs.sq)/n
    (n/2)*log(sigma.sq) + (n/2)*log(2*pi) + (1/(2*sigma.sq))*sum(diffs.sq)
  }
  
  # Calculate metabolism by non linear minimization of an MLE function
  river.mle <- with(
    data, 
    nlm(onestation_negloglik, p=c(GPP=3,ER=-5), 
        DO.obs=DO.obs, DO.deficit=DO.deficit, depth=depth, k.O2=k.O2, 
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
  get_data(metab_model) %>%
    mutate(date=as.Date(format(date.time, "%Y-%m-%d"))) %>%
    group_by(date) %>%
    do(with(., {
      # get today's metabolism estimates
      GPP <- metab_ests[metab_ests$date==date[1], "GPP"]
      ER <- metab_ests[metab_ests$date==date[1], "ER"]

      # prepare other data
      timestep.days <- unique(as.numeric(diff(date.time), units="days"))
      frac.GPP=light/sum(light)
      frac.ER=timestep.days
      frac.D=timestep.days
      
      # model DO
      DO.mod <- numeric(length(DO.obs))
      DO.mod[1] <- DO.obs[1]
      DO.delta <- 
        GPP * frac.GPP / depth + 
        ER * frac.ER / depth + 
        k.O2 * DO.deficit * frac.D
      DO.mod[-1] <- DO.mod[1] + cumsum(DO.delta[-1])
      
      data.frame(., DO.mod=DO.mod)
    }))
    
}
