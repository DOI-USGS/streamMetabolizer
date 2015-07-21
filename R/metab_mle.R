#' @include metab_model-class.R
NULL

#' Maximum likelihood metabolism model fitting function
#' 
#' Fits a model to estimate GPP and ER from input data on DO, temperature, 
#' light, etc.
#' 
#' @param data data.frame with columns having the same names, units, and format 
#'   as the default. See \code{\link{mm_data}} for a full data description.
#' @param info Any metadata you would like to package within the metabolism 
#'   model.
#' @inheritParams mm_is_valid_day
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER, 
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @return A metab_mle object containing the fitted model.
#'   
#' @import dplyr
#' @author Alison Appling, Jordan Read; modeled on LakeMetabolizer
#' @examples
#' \dontrun{
#'  metab_mle(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_mle <- function(
  data=mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light),
  calc_DO_fun=calc_DO_mod) {
  info=NULL, # args for new('metab_mle')
  tests=c('full_day', 'even_timesteps', 'complete_data'), day_start=-1.5, day_end=30 # args for mm_is_valid_day, mm_model_by_ply
) {
  
  # Check data for correct column names & units
  data <- mm_validate_data(data, "metab_mle")
  
  # model the data, splitting into overlapping 31.5-hr 'plys' for each date
  mle.all <- mm_model_by_ply(
    data, mle_1ply, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
    tests=tests, # for mm_is_valid_day
  
  # Package and return results
  new("metab_mle", 
      info=info,
      fit=mle.all,
      args=list(K600=K600, calc_DO_fun=calc_DO_fun, tests=tests, day_start=day_start, day_end=day_end), # keep in order passed to function
      data=data,
      pkg_version=as.character(packageVersion("streamMetabolizer")))
}


#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_mle().
#' 
#' @param data_ply data.frame of the form \code{mm_data(local.time, DO.obs, DO.sat,
#'   depth, temp.water, light)} and containing data for just one estimation-day
#'   (this may be >24 hours but only yields estimates for one 24-hour period)
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER,
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
#' @keywords internal
mle_1ply <- function(data_ply, calc_DO_fun=calc_DO_mod) {
                     tests=tests, day_start=day_start, day_end=day_end, ...) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  timestep.days <- suppressWarnings(mean(as.numeric(diff(v(data_ply$local.time)), units="days"), na.rm=TRUE))
  stop_strs <- mm_is_valid_day(
    data_ply, # data split by mm_model_by_ply
    tests=tests, day_start=day_start, day_end=day_end, # args passed from metab_mle
    timestep_days=timestep.days, need_complete=c("DO.obs","DO.sat","depth","temp.water","light")) # args supplied here
  warn_strs <- character(0)

  # Calculate metabolism by non linear minimization of an MLE function
  if(length(stop_strs) == 0) {
    timestep.days <- suppressWarnings(mean(as.numeric(diff(v(data_ply$local.time)), units="days"), na.rm=TRUE))
    date <- names(which.max(table(as.Date(data_ply$local.time))))
    frac.GPP <- data_ply$light/sum(data_ply[strftime(data_ply$local.time,"%Y-%m-%d")==date,'light'])
    mle.1d <- tryCatch({
      # first: try to run the MLE fitting function
      nlm(onestation_negloglik, p=c(GPP=3, ER=-5, K600=5), 
          DO.obs=data_ply$DO.obs, DO.sat=data_ply$DO.sat, depth=data_ply$depth, temp.water=data_ply$temp.water,
          frac.GPP=frac.GPP, frac.ER=timestep.days, frac.D=timestep.days,
          calc_DO_fun=calc_DO_fun)
    }, warning=function(war) {
      # on warning: record the warning and run nlm again
      warn_strs <- c(warn_strs, war$message)
      suppressWarnings(
        nlm(onestation_negloglik, p=c(GPP=3, ER=-5, K600=5), 
            DO.obs=data_ply$DO.obs, DO.sat=data_ply$DO.sat, depth=data_ply$depth, temp.water=data_ply$temp.water,
            frac.GPP=frac.GPP, frac.ER=timestep.days, frac.D=timestep.days,
            calc_DO_fun=calc_DO_fun)
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
#' Called from mle_1ply(). From ?nlm, this function should be "the function
#' to be minimized, returning a single numeric value. This should be a function 
#' with first argument a vector of the length of p followed by any other 
#' arguments specified by the ... argument."
#' 
#' @param params a vector of length 2, where the first element is GPP and the 
#'   second element is ER (both mg/L/d)
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER,
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @keywords internal
negloglik_1ply <- function(params, K600.daily, DO.obs, DO.sat, depth, temp.water, 
                           frac.GPP, frac.ER, frac.D, calc_DO_fun) {
  
  # Count how many DO observations/predictions there should be
  n <- length(DO.obs)
  
  # Parse params vector (passed from nlm) and produce DO.mod estimates. It's
  # important that these arguments are named because various calc_DO_funs take 
  # slightly different subsets of these arguments
  DO.mod <- calc_DO_fun(
    GPP.daily=params[1], ER.daily=params[2], K600.daily=params[3],
    DO.obs=DO.obs, DO.sat=DO.sat, depth=depth, temp.water=temp.water, 
    frac.GPP=frac.GPP, frac.ER=frac.ER, frac.D=frac.D, DO.mod.1=DO.obs[1], n=n)
  
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
#' @inheritParams predict_DO
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @export
#' @family predict_DO
predict_DO.metab_mle <- function(metab_model) {
  
  # pull args from the model
  calc_DO_fun <- get_args(metab_model)$calc_DO_fun
  day_start <- get_args(metab_model)$day_start
  day_end <- get_args(metab_model)$day_end
  
  # get the metabolism (GPP, ER) data and estimates
  metab_ests <- predict_metab(metab_model)
  data <- get_data(metab_model)
  
  # re-process the input data with the metabolism estimates to predict DO
  mm_model_by_ply(
    data=data, model_fun=mm_predict_1ply, day_start=day_start, day_end=day_end, # for mm_model_by_ply
    calc_DO_fun=calc_DO_fun, metab_ests=metab_ests) # for mm_predict_1ply

}
