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
#'   \item{ \code{temp.water} water temperature, \eqn{degC}}.
#'   
#'   \item{ \code{light} photosynthetically active radiation, \eqn{\mu mol\ 
#'   m^{-2} s^{-1}}{micro mols / m^2 / s}}
#'   
#'   }
#' @param ... additional arguments
#' @return A metab_simple object containing the fitted model.
#' 
#' @import dplyr  
#' @author Alison Appling, Jordan Read; modeled on lakeMetabolizer
#' @examples
#' \dontrun{
#'  metab_simple(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_simple <- function(data, ...) {
  
  # Check data for correct column names
  expected.colnames <- c("date.time","DO.obs","DO.sat","depth","temp.water","light")
  if(!all(expected.colnames %in% names(data))) {
    stop(paste0("data must contain (at least) columns with the names ", paste0(expected.colnames, collapse=", ")))
  }
  
  #' Return the likelihood value for a given set of parameters and observations
  #' 
  #' From ?nlm, this function should be "the function to be minimized, returning
  #' a single numeric value. This should be a function with first argument a 
  #' vector of the length of p followed by any other arguments specified by the 
  #' ... argument."
  #' 
  #' @param params a vector of length 2, where the first element is GPP and the second element is ER (both mg/L/d)
  onestation_negloglik <- function(params, DO.obs, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D) {

    # Count how many DO observations/predictions there should be
    n <- length(DO.obs)
    
    # Parse params vector (passed from nlm) and produce DO.mod estimates
    DO.mod <- .core_model_metab_simple(
      GPP.daily=params[1], ER.daily=params[2], k.600.daily=params[3],
      DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.obs[1], n)
    
    # calculate & return the negative log likelihood of DO.mod values relative
    # to DO.obs values. equivalent to Bob's original code & formula at
    # http://www.statlect.com/normal_distribution_maximum_likelihood.htm
    diffs.sq <- (DO.obs-DO.mod)^2 
    sigma.sq <- sum(diffs.sq)/n
    (n/2)*log(sigma.sq) + (n/2)*log(2*pi) + (1/(2*sigma.sq))*sum(diffs.sq)
  }
  
  # define function to make daily metabolism estimates from the input data
  mle.dummy <- list(minimum=NA, estimate=c(NA,NA,NA), gradient=c(NA,NA,NA), code=NA, iterations=NA)
  est_metab_1d <- function(day) {
    
    # Provide ability to skip a poorly-formatted day for calculating 
    # metabolism, without breaking the whole loop. Just collect 
    # problems/errors as a list of strings and proceed. Also collect warnings.
    stop_strs <- warn_strs <- list()
    
    ## Error checks:
    # Require that the data consist of three consecutive days (10:30 pm on Day 1 to 6 am on Day 3)
    if(!isTRUE(all.equal(diff(range(day$date.time)) %>% as.numeric(units="days"), 
                         as.difftime(31.5, units="hours") %>% as.numeric(units="days"), 
                         tol=as.difftime(31, units="mins") %>% as.numeric(units="days")))) {
      stop_strs <- c(stop_strs, "incomplete time series")
    }
    # Require that on each day date.time has a ~single, ~consistent time step
    timestep.days <- suppressWarnings(mean(as.numeric(diff(day$date.time), units="days"), na.rm=TRUE))
    timestep.deviations <- suppressWarnings(diff(range(as.numeric(diff(day$date.time), units="days"), na.rm=TRUE)))
    if(length(stop_strs) == 0 & is.finite(timestep.days) & is.finite(timestep.deviations)) {
      # max-min timestep length can't be more than 1% of mean timestep length
      if((timestep.deviations / timestep.days) > 0.001) { 
        stop_strs <- c(stop_strs, "uneven timesteps")
      }
      # all timesteps per day must add up to 31.5 hrs (31.5/24 days), plus or minus 0.51 hrs
      if(abs(timestep.days * length(day$date.time) - 31.5/24) > 0.51/24) { 
        stop_strs <- c(stop_strs, paste0("sum(timesteps) != 31.5 hours"))
      }
    } else {
      stop_strs <- c(stop_strs, "can't measure timesteps")
    }
    # Require complete data
    if(any(is.na(day$DO.obs))) stop_strs <- c(stop_strs, "NAs in DO.obs")
    if(any(is.na(day$DO.sat))) stop_strs <- c(stop_strs, "NAs in DO.sat")
    if(any(is.na(day$depth))) stop_strs <- c(stop_strs, "NAs in depth")
    if(any(is.na(day$temp.water))) stop_strs <- c(stop_strs, "NAs in temp.water")
    if(any(is.na(day$light))) stop_strs <- c(stop_strs, "NAs in light")
    
    # Calculate metabolism by non linear minimization of an MLE function
    if(length(stop_strs) == 0) {
      mle.1d <- tryCatch({
        # first: try to run the MLE fitting function
        nlm(onestation_negloglik, p=c(GPP=3, ER=-5, k.600=5), 
            DO.obs=day$DO.obs, DO.sat=day$DO.sat, depth=day$depth, temp.water=day$temp.water,
            frac.GPP=day$light/sum(day$light), frac.ER=timestep.days, frac.D=timestep.days)
      }, warning=function(war) {
        # on warning: record the warning and run nlm again
        warn_strs <- c(warn_strs, war$message)
        suppressWarnings(
          nlm(onestation_negloglik, p=c(GPP=3, ER=-5, k.600=5), 
              DO.obs=day$DO.obs, DO.sat=day$DO.sat, depth=day$depth, temp.water=day$temp.water,
              frac.GPP=day$light/sum(day$light), frac.ER=timestep.days, frac.D=timestep.days)
        )
      }, error=function(err) {
        # on error: give up, remembering error. dummy values provided below
        stop_strs <- c(stop_strs, err$message)
      })
    } 
    
    # Check again - stop_strs may have accumulated during nlm() call. If
    # stopped, record why.
    if(length(stop_strs) > 0) {
      mle.1d <- mle.dummy
    }
    
    # Return
    data.frame(GPP=mle.1d$estimate[1], ER=mle.1d$estimate[2], K600=mle.1d$estimate[3],
               grad.GPP=mle.1d$gradient[1], grad.ER=mle.1d$gradient[2], grad.K600=mle.1d$gradient[3],
               minimum=mle.1d$minimum, code=mle.1d$code, iterations=mle.1d$iterations, 
               warnings=paste0(unlist(warn_strs), collapse="; "), 
               errors=paste0(unlist(stop_strs), collapse="; "),
               stringsAsFactors=FALSE)
  }
  
  # Identify the data plys that will let us use a 31.5-hr window for each date -
  # this labeling can be stored in two additional columns (odd.- and even.- 
  # date.group)
  data.plys <- data %>% 
    mutate(date=as.Date(format(date.time, "%Y-%m-%d")),
           hour=24*(convert_date_to_doyhr(date.time) %% 1))
  unique.dates <- unique(data.plys$date)
  odd.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 1)]
  even.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 0)]
  data.plys <- data.plys %>% 
    group_by(date) %>%
    mutate(odd.date.group=if(date[1] %in% odd.dates) date else c(date[1]-1, as.Date(NA), date[1]+1)[ifelse(hour <= 6, 1, ifelse(hour < 22.5, 2, 3))],
           even.date.group=if(date[1] %in% even.dates) date else c(date[1]-1, as.Date(NA), date[1]+1)[ifelse(hour <= 6, 1, ifelse(hour < 22.5, 2, 3))]) %>%
    ungroup() %>% select(-date)
  
  # Make daily metabolism estimates for each ply of the data, using two
  # group_by/do combinations to cover the odd and even groupings
  . <- ".dplyr.var"
  mle.all <- 
    bind_rows(
      data.plys %>% group_by(date=odd.date.group) %>% do(est_metab_1d(.)), # filter(!is.na(odd.date.group)) %>%, filter(!is.na(even.date.group)) %>% 
      data.plys %>% group_by(date=even.date.group) %>% do(est_metab_1d(.))) %>% 
    filter(!is.na(date), date %in% unique.dates) %>%
    arrange(date) 
  
  # Package and return results. There are no args, just data, for this
  # particular model
  new("metab_simple", 
      fit=mle.all,
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
.core_model_metab_simple <- function(GPP.daily, ER.daily, k.600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.obs.1, n) {
  # partition GPP and ER into their timestep-specific rates (mg/L/timestep at
  # each timestep)
  GPP <- GPP.daily * frac.GPP / depth
  ER <- ER.daily * frac.ER / depth
  K <- convert_k600_to_kGAS(k.600.daily, temperature=temp.water, gas="O2") * frac.D
  
  # make sure anything in the following loop has n observations
  if(length(DO.sat) != n) DO.sat <- rep(DO.sat, length.out=n) # this would be odd
  
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.obs.1
  for(i in 2:n) {
    DO.mod[i] <- DO.mod[i-1] +
      GPP[i] + 
      ER[i] + 
      K[i] * (DO.sat[i] - DO.mod[i-1])
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
#' @import dplyr
#' @family predict_metab
predict_metab.metab_simple <- function(metab_model) {
  
  get_fit(metab_model) %>% 
    mutate(NEP=GPP+ER) %>%
    select(date, GPP, ER, NEP, K600)
  
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
      # get the daily metabolism estimates, and skip today if they're missing
      metab_est <- metab_ests[metab_ests$date==date[1],]
      if(is.na(metab_est$GPP)) {
        return(data.frame(., DO.mod=NA))
      }
      
      # prepare auxiliary data
      n <- length(date.time)
      timestep.days <- mean(as.numeric(diff(date.time), units="days"))
      frac.GPP=light/sum(light)
      frac.ER=timestep.days
      frac.D=timestep.days
      
      # produce DO.mod estimates for today's GPP and ER
      DO.mod <- .core_model_metab_simple(
        GPP.daily=metab_est$GPP, 
        ER.daily=metab_est$ER, 
        k.600.daily=metab_est$K600, 
        DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.obs[1], n)
      
      data.frame(., DO.mod=DO.mod)
    }))
  
}
