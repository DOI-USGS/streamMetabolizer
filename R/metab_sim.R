#' @include metab_model-class.R
NULL

#' Simulate dissolved oxygen data from input data
#' 
#' Takes input data
#' 
#' @author Alison Appling, Maite Arroita
#' @inheritParams metab_model_prototype
#' @inheritParams mm_is_valid_day
#' @return A metab_sim object containing the fitted model.
#' @examples
#' \dontrun{
#'  metab_sim(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_sim <- function(
  data=mm_data(local.time, DO.sat, depth, temp.water, light), data_daily=mm_data(local.date, GPP, ER, K600), info=NULL, day_start=-1.5, day_end=30, # inheritParams metab_model_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # args for mm_is_valid_day
  model_specs=list(calc_DO_fun=calc_DO_mod)
) {
  
  # Check data for correct column names & units
  dat_list <- mm_validate_data(data, if(missing(data_daily)) NULL else data_daily, "metab_sim")
  data <- dat_list[['data']]
  data_daily <- dat_list[['data_daily']]
  
  # there could be a model_by_ply call in here to copy data for overlapping
  # days, but maybe we can get away with doing nothing until the predict phase.
  sim_all <- data_daily # GPP, ER, etc. were given as data but will become our predictors
  
  # Package and return results
  metab_model(
    model_class="metab_sim", 
    info=info,
    fit=sim_all,
    args=list(day_start=day_start, day_end=day_end, tests=tests),
    data=data,
    data_daily=data_daily)
}


#### helpers ####



#### metab_sim class ####

#' Data simulator
#' 
#' \code{metab_sim} models generate a DO time series from other input data,
#' including GPP, ER, and K600 values
#' 
#' @exportClass metab_sim
#' @family metab.model.classes
setClass(
  "metab_sim", 
  contains="metab_model"
)


#' Explain why we can't make dissolved oxygen predictions from a metab_sim.
#' 
#' metab_sim only fits ER and K, and only for the darkness hours. While it 
#' would be possible to make predictions just for those hours, it'd be costly to
#' implement and has not-yet-obvious benefits.
#' 
#' @inheritParams predict_DO
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @export
#' @family predict_DO
predict_DO.metab_sim <- function(metab_model, ...) {
  
  # pull args from the model
  calc_DO_fun <- calc_DO_mod # this isn't used in the model fitting but makes sense for prediction
  day_start <- get_args(metab_model)$day_start
  day_end <- get_args(metab_model)$day_end
  
  # get the metabolism (GPP, ER) data and estimates
  metab_ests <- predict_metab(metab_model)
  data <- get_data(metab_model)
  
  # re-process the input data with the metabolism estimates to predict DO, using
  # our special nighttime regression prediction function
  mm_model_by_ply(
    model_fun=metab_sim_predict_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, # for mm_model_by_ply
    calc_DO_fun=calc_DO_fun) # for mm_predict_1ply
  
}

#' Helper to predict_DO.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER, 
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @return a data.frame of predictions
#' @importFrom stats complete.cases
metab_sim_predict_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, local_date, # inheritParams mm_model_by_ply_prototype
  calc_DO_fun
) {
  warning("unimplemented")
}