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
#' @author Alison Appling, Maite Arroita
#' @inheritParams metab_model_prototype
#' @inheritParams mm_is_valid_day
#' @return A metab_sim object containing the fitted model.
#' @examples
#' # set the date in several formats
#' start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
#' end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
#' start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="Etc/GMT+7")
#' end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="Etc/GMT+7")
#' mid.date <- as.Date(start.posix + (end.posix - start.posix)/2)
#' start.numeric <- as.numeric(start.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'    tz="Etc/GMT+7"), units='hours')
#' end.numeric <- as.numeric(end.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'   tz="Etc/GMT+7"), units='hours')
#' 
#' # get, format, & subset data
#' vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
#' vfrenchshort <- vfrench[vfrench$local.time >= start.posix & vfrench$local.time <= end.posix, ]
#' vdaily <- data.frame(local.date="2012-08-24", GPP=2, ER=-3, K600=21, stringsAsFactors=FALSE)
#' 
#' # sim
#' mm <- metab_sim(data=vfrenchshort, data_daily=vdaily,
#'   day_start=start.numeric, day_end=end.numeric, 
#'   model_specs=specs('s_np_oipcpi_eu_.rnorm', err.proc.sigma=0.07))
#' get_fit(mm)
#' get_data_daily(mm)
#' get_fitting_time(mm)
#' plot_DO_preds(predict_DO(mm))
#' plot_DO_preds(predict_DO(mm))
#' \dontrun{
#'  metab_sim(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_sim <- function(
  data=mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light, optional='DO.obs'), # inheritParams metab_model_prototype
  data_daily=mm_data(local.date, DO.mod.1, GPP, ER, K600, optional='DO.mod.1'), # inheritParams metab_model_prototype
  model_specs=specs('s_np_oipcpi_eu_.rnorm'), # inheritParams metab_model_prototype
  info=NULL, day_start=4, day_end=27.99, # inheritParams metab_model_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data') # args for mm_is_valid_day
  
) {
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(data, if(missing(data_daily)) NULL else data_daily, "metab_sim")
    
    # Move the simulation-relevant parameters to calc_DO_args for use in predict_DO
    calc_DO_arg_names <- c('err.obs.sigma','err.obs.phi','err.proc.sigma','err.proc.phi','ODE_method')
    model_specs$calc_DO_args = model_specs[calc_DO_arg_names]
    model_specs <- model_specs[-which(names(model_specs) %in% calc_DO_arg_names)]
  })
  
  # Package and return results
  metab_model(
    model_class="metab_sim", 
    info=info,
    fit=dat_list[['data_daily']], # GPP, ER, etc. were given as data but will become our predictors
    fitting_time=fitting_time,
    args=list(model_specs=model_specs, day_start=day_start, day_end=day_end, tests=tests),
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
#' Makes daily predictions of GPP, ER, and NEP.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @importFrom stats qnorm setNames
#' @export
#' @family predict_metab
predict_metab.metab_sim <- function(metab_model, ...) {
  
  # Select only those columns required for metabolism prediction and available
  # in the fit. At present this is all of the columns, but that could change
  fit <- get_fit(metab_model)
  vars <- c("local.date","DO.mod.1","GPP","ER","K600")
  fit[vars[vars %in% names(fit)]]
  
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
predict_DO.metab_sim <- function(metab_model, ...) {

  # call the generic, which generally does what we want
  sim.seed <- get_args(metab_model)$model_specs$sim.seed
  if(!is.na(sim.seed)) set.seed(sim.seed)
  preds_w_err <- NextMethod(calc_DO_fun=calc_DO_mod_w_sim_error)
  preds_wo_err <- NextMethod(calc_DO_fun=calc_DO_mod)
  
  # copy the predictions into the DO.obs column to complete the simulation
  preds_wo_err$DO.obs <- preds_w_err$DO.mod
  
  preds_wo_err
}
