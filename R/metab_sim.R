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
#' @importFrom unitted v
#' @examples
#' # start with non-DO data (DO used only to pick first DO of each day)
#' dat <- data_metab('3', res='15')
#' dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)),
#'   GPP.daily=2, ER.daily=-3, K600.daily=21, stringsAsFactors=FALSE)
#' 
#' # define simulation parameters
#' mm <- metab_sim(
#'   specs(mm_name('sim'), err_obs_sigma=0.1, err_proc_sigma=2),
#'   data=dat, data_daily=dat_daily)
#' # actual simulation happens during prediction - different each time
#' predict_DO(mm)[seq(1,50,by=10),]
#' predict_DO(mm)[seq(1,50,by=10),]
#' 
#' # or same each time if seed is set
#' mm <- metab_sim(
#'   specs(mm_name('sim'), err_obs_sigma=0.1, err_proc_sigma=2, sim_seed=248),
#'   data=dat, data_daily=dat_daily)
#' predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#' predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#' 
#' # fancy GPP equation
#' dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)),
#'   Pmax=8, alpha=0.01, ER.daily=-3, K600.daily=21, stringsAsFactors=FALSE)
#' mm <- metab_sim(
#'   specs(mm_name('sim', GPP_fun='satlight'), err_obs_sigma=0.1, err_proc_sigma=2),
#'   data=dat, data_daily=dat_daily)
#' get_params(mm)
#' predict_metab(mm) # metab estimates are for data without errors
#' predict_DO(mm)[seq(1,50,by=10),]
#' 
#' \dontrun{
#' plot_DO_preds(predict_DO(mm))
#' plot_DO_preds(mm)
#' plot_DO_preds(mm, y_var='conc') + geom_line(aes(y=DO.pure), color='tan', alpha=0.8, size=1)
#' }
#' @export
#' @family metab_model
metab_sim <- function(
  specs=specs(mm_name('sim')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light, optional='DO.obs'),
  data_daily=mm_data(date, DO.mod.1, GPP.daily, Pmax, alpha, ER.daily, ER20, K600.daily, optional='all'),
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
  })
  
  # Package and return results
  metab_model(
    model_class="metab_sim", 
    info=info,
    fit=NULL,
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
#' @importFrom stats rnorm
#' @family predict_DO
predict_DO.metab_sim <- function(metab_model, date_start=NA, date_end=NA, ...) {

  # simulate errors to add to modeled data
  specs <- get_specs(metab_model)
  sim_seed <- specs$sim_seed
  if(!is.na(sim_seed)) set.seed(sim_seed)
  n <- nrow(get_data(metab_model))
  err.obs <- as.numeric(stats::filter(rnorm(n, 0, specs$err_obs_sigma), filter=specs$err_obs_phi, method="recursive"))
  err.proc <- as.numeric(stats::filter(rnorm(n, 0, specs$err_proc_sigma), filter=specs$err_proc_phi, method="recursive"))
  
  # call the generic a few times to get DO with proc and proc+obs error
  preds <- NextMethod(use_saved=FALSE)
  preds$DO.pure <- preds$DO.mod # DO.mod has the DO implied by the daily metab params (error-free, not predicted by anybody)
  metab_model@data$err.proc <- err.proc
  preds$DO.mod <- NextMethod(use_saved=FALSE)$DO.mod # DO.mod has the 'true' DO (with proc err)
  metab_model@data$err.obs <- err.obs
  preds$DO.obs <- NextMethod(use_saved=FALSE)$DO.mod # the 'observed' DO (with obs err)
  
  # add additional observation error in the form of DO rounding if requested
  if(!is.na(specs$err_round)) preds_w_err$DO.obs <- round(preds_w_err$DO.obs, digits=specs$err_round)
  
  # return
  preds
}
