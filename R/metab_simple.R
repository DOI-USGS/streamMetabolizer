#' @title Basic metabolism model fitting function
#' @description Fits a model to estimate GPP and ER from input data on DO, 
#'   temperature, light, etc.
#'   
#' @param data data.frame with columns \itemize{ \item{ DO.obs Vector of 
#'   dissolved oxygen concentration observations, \eqn{mg O[2] L^{-1}}{mg O2 / 
#'   L}} \item{ DO.sat Vector of dissolved oxygen saturation values based on 
#'   water temperature. Calculate using \link{o2.at.sat}} \item{ k.gas Vector of
#'   kGAS values calculated from any of the gas flux models (e.g., 
#'   \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}} \item{ PAR 
#'   Vector of photosynthetically active radiation in \eqn{\mu mol\ m^{-2} 
#'   s^{-1}}{micro mols / m^2 / s}} \item{ temp.water Vector of water 
#'   temperatures in \eqn{^{\circ}C}{degrees C}. Used in scaling respiration 
#'   with temperature} }
#' @param ... additional arguments
#' @return A data.frame with columns corresponding to components of metabolism 
#'   \describe{ \item{GPP}{numeric estimate of Gross Primary Production, \eqn{mg
#'   O_2 L^{-1} d^{-1}}{mg O2 / L / d}} \item{R}{numeric estimate of
#'   Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}} \item{NEP}{numeric
#'   estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L
#'   / d}} }
#'   
#' @details The model has inputs and parameters
#'   
#' @author Alison Appling, Jordan Read; modeled on lakeMetabolizer
#' @examples
#' \dontrun{
#'  metab_simple(data=data.frame(empty="shouldbreak"))
#'  
#'  # use a subset of data from Bob
#'  library(dplyr)
#'  french_data <- read.csv("data/french_data.csv") %>% mutate(date.time=as.POSIXct(strptime(date.time, format="%Y-%m-%d %H:%M:%S")))
#'  french_units <- as.character(c("", read.csv("data/french_units.csv")$x))
#'  library(unitted); french <- u(french_data, french_units)
#'  mm <- metab_simple(data=v(french))
#'  slot(mm, "fit")
#' }
#' @export
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
    
    # parse params vector (passed from nlm)
    GPP <- params[1]
    ER <- params[2]
    
    # model DO with given params
    DO.mod <- numeric(length(DO.obs))
    DO.mod[1] <- DO.obs[1]
    DO.delta <- 
      GPP * frac.GPP / depth + 
      ER * frac.ER / depth + 
      k.O2 * DO.deficit * frac.D # Must we use DO.mod[t-1] instead of DO.obs[t-1] (i.e., from DO.deficit) here?
    DO.mod[-1] <- DO.mod[1]+cumsum(DO.delta[-1])
    
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
  
  metab_model(fit=river.mle)
}