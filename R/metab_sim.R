#' @include metab_model-class.R
NULL

#' Simulate dissolved oxygen data from input data
#' 
#' Takes input data in the form of a sub-daily time series (\code{data}) of 
#' DO.sat, depth, temperature, and light, and a daily time series 
#' (\code{data_daily}) of GPP, ER, and K600 values, and turns these into 
#' simulated DO.obs. Either \code{data} or \code{data_daily} should specify a 
#' starting DO.obs value for each day; if in \code{data}, this takes the form of
#' a DO.obs column with values on at least the first time point of each day (all
#' other values are ignored), or if in \code{data_daily}, this takes the form of
#' a DO.mod.1 column with one starting DO value per day.
#' 
#' @author Alison Appling, Bob Hall
#'   
#' @inheritParams metab
#' @return A metab_sim object containing the fitted model. This object can be 
#'   inspected with the functions in the \code{\link{metab_model_interface}}.
#'   
#' @examples
#' # start with non-DO data (DO used only to pick first DO of each day)
#' dat <- data_metab('3', res='15')
#' dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)),
#'   GPP=2, ER=-3, K600=21, stringsAsFactors=FALSE)
#' 
#' # define simulation parameters
#' mm <- metab_sim(
#'   specs(mm_name('sim'), err.obs.sigma=0.01, err.proc.sigma=0.1),
#'   data=dat, data_daily=dat_daily)
#' # actual simulation happens during prediction - different each time
#' predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#' predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#' 
#' # or same each time if seed is set
#' mm <- metab_sim(
#'   specs(mm_name('sim'), err.obs.sigma=0.01, err.proc.sigma=0.1, sim.seed=248),
#'   data=dat, data_daily=dat_daily)
#' predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#' predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#' 
#' \dontrun{
#' plot_DO_preds(predict_DO(mm))
#' plot_DO_preds(predict_DO(mm))
#' }
#' @export
#' @family metab_model
metab_sim <- function(
  specs=specs(mm_name('sim')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light, optional='DO.obs'),
  data_daily=mm_data(date, DO.mod.1, GPP, ER, K600, optional='DO.mod.1'),
  info=NULL
) {
  
  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name('sim'))
  }
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, if(missing(data_daily)) NULL else data_daily, "metab_sim")
    
    # Move the simulation-relevant parameters to calc_DO_args for use in predict_DO
    calc_DO_arg_names <- c('err.obs.sigma','err.obs.phi','err.proc.sigma','err.proc.phi','ODE_method')
    specs$calc_DO_args = specs[calc_DO_arg_names]
    specs <- specs[-which(names(specs) %in% calc_DO_arg_names)]
  })
  
  # Package and return results
  metab_model(
    model_class="metab_sim", 
    info=info,
    fit=dat_list[['data_daily']], # GPP, ER, etc. were given as data but will become our predictors
    fitting_time=fitting_time,
    specs=specs,
    data=dat_list[['data']],
    data_daily=dat_list[['data_daily']])
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


#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily predictions of GPP, ER, and K600.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @import dplyr
#' @export
#' @family predict_metab
predict_metab.metab_sim <- function(metab_model, date_start=NA, date_end=NA, ...) {
  
  # Select only those columns required for metabolism prediction and available
  # in the fit. At present this is all of the columns, but that could change
  fit <- get_fit(metab_model) %>%
    mm_filter_dates(date_start=date_start, date_end=date_end)
  vars <- c("date","DO.mod.1","GPP","ER","K600")
  fit[vars[vars %in% names(fit)]] %>%
    mutate(warnings=as.character(NA), errors=as.character(NA))
}


#' Simulate values for DO.obs (with error) and DO.mod (without)
#' 
#' Generate simulated values for DO.obs (including any error specified in the 
#' call to \code{metab_sim()}) and DO.mod (with no error, just the predictions 
#' from the specified GPP, ER, and K600). The errors are randomly generated on 
#' every new call to predict_DO.
#' 
#' @inheritParams predict_DO
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @export
#' @family predict_DO
predict_DO.metab_sim <- function(metab_model, date_start=NA, date_end=NA, ...) {

  # call the generic, which generally does what we want
  sim.seed <- get_specs(metab_model)$sim.seed
  if(!is.na(sim.seed)) set.seed(sim.seed)
  preds_w_err <- NextMethod(calc_DO_fun=calc_DO_mod_w_sim_error, use_saved=FALSE)
  preds_wo_err <- NextMethod(calc_DO_fun=calc_DO_mod, use_saved=FALSE)
  
  # copy the predictions into the DO.obs column to complete the simulation
  preds_wo_err$DO.obs <- preds_w_err$DO.mod
  
  preds_wo_err
}
