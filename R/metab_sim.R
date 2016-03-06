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
#' # set the date in several formats
#' start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
#' end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
#' start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' mid.date <- as.Date(start.posix + (end.posix - start.posix)/2, tz=lubridate::tz(start.posix))
#' start.numeric <- as.numeric(start.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'    tz="UTC"), units='hours')
#' end.numeric <- as.numeric(end.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'   tz="UTC"), units='hours')
#' 
#' # get, format, & subset data
#' vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
#' vfrenchshort <- vfrench[vfrench$solar.time >= start.posix & vfrench$solar.time <= end.posix, ]
#' vdaily <- data.frame(date="2012-08-24", GPP=2, ER=-3, K600=21, stringsAsFactors=FALSE)
#' 
#' # sim
#' mm <- metab_sim(specs=specs('s_np_oipcpi_eu_.rnorm', err.proc.sigma=0.07, 
#'   day_start=start.numeric, day_end=end.numeric), data=vfrenchshort, data_daily=vdaily)
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
  specs=specs(mm_name('sim')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light, optional='DO.obs'),
  data_daily=mm_data(date, DO.mod.1, GPP, ER, K600, optional='DO.mod.1'),
  info=NULL
) {
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(data, if(missing(data_daily)) NULL else data_daily, "metab_sim")
    
    # Move the simulation-relevant parameters to calc_DO_args for use in predict_DO
    calc_DO_arg_names <- c('err.obs.sigma','err.obs.phi','err.proc.sigma','err.proc.phi','ODE_method','day_start','day_end','tests')
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
#' @export
#' @family predict_metab
predict_metab.metab_sim <- function(metab_model, date_start=NA, date_end=NA, ...) {
  
  # Select only those columns required for metabolism prediction and available
  # in the fit. At present this is all of the columns, but that could change
  fit <- get_fit(metab_model) %>%
    mm_filter_dates(date_start=date_start, date_end=date_end)
  vars <- c("date","DO.mod.1","GPP","ER","K600")
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
predict_DO.metab_sim <- function(metab_model, date_start=NA, date_end=NA, ...) {

  # call the generic, which generally does what we want
  sim.seed <- get_specs(metab_model)$sim.seed
  if(!is.na(sim.seed)) set.seed(sim.seed)
  preds_w_err <- NextMethod(calc_DO_fun=calc_DO_mod_w_sim_error)
  preds_wo_err <- NextMethod(calc_DO_fun=calc_DO_mod)
  
  # copy the predictions into the DO.obs column to complete the simulation
  preds_wo_err$DO.obs <- preds_w_err$DO.mod
  
  preds_wo_err
}
