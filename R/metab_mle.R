#' @include metab_model-class.R
NULL

#' Basic metabolism model fitting function
#' 
#' Fits a model to estimate GPP and ER from input data on DO, 
#'   temperature, light, etc.
#'   
#' @param data data.frame with columns having the same names, units, and format 
#'   as the default. See \code{\link{mm_data}} for a full data description.
#' @param ... additional arguments
#' @return A metab_mle object containing the fitted model.
#'   
#' @import dplyr
#' @author Alison Appling, Jordan Read; modeled on lakeMetabolizer
#' @examples
#' \dontrun{
#'  metab_mle(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_mle <- function(
  data=mm_data(date.time, DO.obs, DO.sat, depth, temp.water, light),
  ...) {
  
  # Check data for correct column names & units
  data <- mm_validate_data(data, "metab_mle")
  
  # Identify the data plys that will let us use a 31.5-hr window for each date -
  # this labeling can be stored in two additional columns (odd.- and even.- 
  # date.group)
  date.time <- hour <- ".dplyr.var"
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
  
  # Estimate daily metabolism for each ply of the data, using two group_by/do
  # combinations to cover the odd and even groupings
  . <- odd.date.group <- even.date.group <- ".dplyr.var"
  mle.all <- 
    bind_rows(
      data.plys %>% group_by(date=odd.date.group) %>% do(est_metab_1d(.)), # filter(!is.na(odd.date.group)) %>%, filter(!is.na(even.date.group)) %>% 
      data.plys %>% group_by(date=even.date.group) %>% do(est_metab_1d(.))) %>% 
    filter(!is.na(date), date %in% unique.dates) %>%
    arrange(date) 
  
  # Package and return results. There are no args, just data, for this 
  # particular model
  new("metab_mle", 
      fit=mle.all,
      args=list(),
      data=data,
      pkg_version=as.character(packageVersion("streamMetabolizer")))
}

#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_mle().
#' 
#' @param day data.frame of the form \code{mm_data(date.time, DO.obs, DO.sat,
#'   depth, temp.water, light)} and containing data for just one estimation-day
#'   (this may be >24 hours but only yields estimates for one 24-hour period)
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
#' @keywords internal
est_metab_1d <- function(day) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  stop_strs <- mm_is_valid_day(day, need_complete=c("DO.obs","DO.sat","depth","temp.water","light"))
  warn_strs <- character(0)

  # Calculate metabolism by non linear minimization of an MLE function
  if(length(stop_strs) == 0) {
    timestep.days <- suppressWarnings(mean(as.numeric(diff(v(day$date.time)), units="days"), na.rm=TRUE))
    mle.1d <- tryCatch({
      # first: try to run the MLE fitting function
      nlm(onestation_negloglik, p=c(GPP=3, ER=-5, K600=5), 
          DO.obs=day$DO.obs, DO.sat=day$DO.sat, depth=day$depth, temp.water=day$temp.water,
          frac.GPP=day$light/sum(day$light), frac.ER=timestep.days, frac.D=timestep.days)
    }, warning=function(war) {
      # on warning: record the warning and run nlm again
      warn_strs <- c(warn_strs, war$message)
      suppressWarnings(
        nlm(onestation_negloglik, p=c(GPP=3, ER=-5, K600=5), 
            DO.obs=day$DO.obs, DO.sat=day$DO.sat, depth=day$depth, temp.water=day$temp.water,
            frac.GPP=day$light/sum(day$light), frac.ER=timestep.days, frac.D=timestep.days)
      )
    }, error=function(err) {
      # on error: give up, remembering error. dummy values provided below
      stop_strs <- c(stop_strs, err$message)
    })
  } 
  
  # stop_strs may have accumulated during nlm() call. If failed, use dummy data 
  # to fill in the model output with NAs.
  if(length(stop_strs) > 0) {
    mle.1d <- list(minimum=NA, estimate=c(NA,NA,NA), gradient=c(NA,NA,NA), code=NA, iterations=NA)
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(GPP=mle.1d$estimate[1], ER=mle.1d$estimate[2], K600=mle.1d$estimate[3],
             grad.GPP=mle.1d$gradient[1], grad.ER=mle.1d$gradient[2], grad.K600=mle.1d$gradient[3],
             minimum=mle.1d$minimum, code=mle.1d$code, iterations=mle.1d$iterations, 
             warnings=paste0(warn_strs, collapse="; "), 
             errors=paste0(stop_strs, collapse="; "),
             stringsAsFactors=FALSE)
}

#' Return the likelihood value for a given set of parameters and observations
#' 
#' Called from est_metab_1d(). From ?nlm, this function should be "the function
#' to be minimized, returning a single numeric value. This should be a function 
#' with first argument a vector of the length of p followed by any other 
#' arguments specified by the ... argument."
#' 
#' @param params a vector of length 2, where the first element is GPP and the 
#'   second element is ER (both mg/L/d)
#' @keywords internal
onestation_negloglik <- function(params, DO.obs, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D) {
  
  # Count how many DO observations/predictions there should be
  n <- length(DO.obs)
  
  # Parse params vector (passed from nlm) and produce DO.mod estimates
  DO.mod <- calc_DO_mod(
    GPP.daily=params[1], ER.daily=params[2], K600.daily=params[3],
    DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.obs[1], n)
  
  # calculate & return the negative log likelihood of DO.mod values relative
  # to DO.obs values. equivalent to Bob's original code & formula at
  # http://www.statlect.com/normal_distribution_maximum_likelihood.htm
  diffs.sq <- (DO.obs-DO.mod)^2 
  sigma.sq <- sum(diffs.sq)/n
  (n/2)*log(sigma.sq) + (n/2)*log(2*pi) + (1/(2*sigma.sq))*sum(diffs.sq)
}





#### metab_mle class ####

#' Metabolism model fitted by maximum likelihood estimation
#' 
#' \code{metab_mle} models use non-linear minimization of the negative log
#' likelihood to fit values of GPP, ER, and K for a given DO curve.
#' 
#' @exportClass metab_mle
#' @family metab.model.classes
setClass(
  "metab_mle", 
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
predict_metab.metab_mle <- function(metab_model) {
  
  GPP <- ER <- NEP <- K600 <- ".dplyr.var"
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
predict_DO.metab_mle <- function(metab_model) {
  
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
      DO.mod <- calc_DO_mod(
        GPP.daily=metab_est$GPP, 
        ER.daily=metab_est$ER, 
        K600.daily=metab_est$K600, 
        DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.obs[1], n)
      
      data.frame(., DO.mod=DO.mod)
    }))
  
}
