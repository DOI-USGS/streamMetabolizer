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
#' ## simulations with variation all at sub-daily scale
#' # prepare input data (DO used only to pick first DO of each day)
#' dat <- data_metab('3', res='15')
#' dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)),
#'   GPP.daily=2, ER.daily=-3, K600.daily=21, stringsAsFactors=FALSE)
#' 
#' # define simulation parameters
#' mm <- metab_sim(
#'   specs(mm_name('sim'), err_obs_sigma=0.1, err_proc_sigma=2),
#'   data=dat, data_daily=dat_daily)
#' # actual simulation happens during prediction - different each time
#' get_params(mm)
#' predict_metab(mm)
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
#' ## simulations with variation at both sub-daily and multi-day scales
#' sp <- specs(mm_name('sim', pool_K600='none'),
#'   K600_daily = 'pmax(0, rnorm(n, 10, 3))') # n is available within sim models
#' mm <- metab(sp, dat, select(dat_daily, -K600.daily))
#' 
#' \dontrun{
#' plot_DO_preds(predict_DO(mm))
#' plot_DO_preds(mm)
#' library(ggplot2)
#' plot_DO_preds(mm, y_var='conc') + geom_line(aes(y=DO.pure), color='tan', alpha=0.8, size=1)
#' }
#' @export
#' @family metab_model
metab_sim <- function(
  specs=specs(mm_name('sim')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light, optional='DO.obs'),
  data_daily=mm_data(date, K600.daily, GPP.daily, Pmax, alpha, ER.daily, ER20, DO.mod.1, optional='all'),
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

#' @describeIn get_params Generates new simulated values for daily parameters if
#'   they were described with evaluatable expressions in \code{\link{specs}}, or
#'   returns the fixed values for daily parameters if they were set in 
#'   \code{data_daily}
#' @importFrom unitted v
#' @export
get_params.metab_sim <- function(
  metab_model, date_start=NA, date_end=NA, 
  uncertainty=c('sd','ci','none'), messages=TRUE, fixed=c('none','columns','stars'), 
  ..., attach.units=FALSE) {
  
  # get model specs, features, and data_daily
  specs <- get_specs(metab_model)
  features <- mm_parse_name(specs$model_name)
  data_daily <-
    if(!is.null(get_data_daily(metab_model))) {
      get_data_daily(metab_model)
    } else {
      # data_daily wasn't given; create a data.frame of dates only based on data
      gimmedates <- function(...) data_frame(ignore=NA) # has to be separate to avoid tripping up deprecation check in convert_date_to_doyhr
      mm_model_by_ply(
        gimmedates, data=get_data(metab_model),
        day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests, timestep_days=FALSE)[1]
    }
  
  # define a controlled environment for evaluation of text to generate
  # parameters
  env <- new.env(parent=emptyenv())
  for(dd in names(data_daily))
  assign(dd, data_daily[dd], envir=env)
  
  # function to get a parameter from data_daily or specs according to
  # availability & format
  get_par <- function(par, required=TRUE) { 
    # option 1: get param from specs, which uses underscore-separated names
    par_ <- gsub("\\.", "_", par)
    if(exists(par_, specs)) {
      parsp <- specs[[par_]]
      if(is.null(parsp)) {
        # option 1a: spec can be NULL as long as there are values in data_daily 
        # instead. keep going to let final check throw an error if it's in
        # neither place
      } else if(is.numeric(parsp[1])) {
        if(length(parsp) == nrow(data_daily)) {
          # option 1b: spec is already numeric; return as-is
          return(parsp)
        } else if(length(parsp) == 1) {
          return(rep(parsp, nrow(data_daily)))
        } else {
          stop(paste0("if numeric, specs$", par_, " must have length 1 or nrow(data_daily)"))
        }
      } else if(is.character(parsp[1]) && length(parsp) == 1) {
        # option 1c: spec is character; return evaluated version. join assumes
        # data_daily and fit have no overlap
        return(eval(parse(text=parsp), envir=env, enclos=emptyenv()))
      } else {
        stop(paste0("specs$", par_, " must be numeric or length-1 character"))
      }
    }
    
    # option 2: if the parameter isn't available in specs, make sure it's 
    # present in data_daily (which uses period-separated names). if it exists,
    # good, and don't add a column to fit (assign NULL instead)
    if(exists(par, data_daily)) {
      return(NULL)
    }
    
    # if we get here then the parameter isn't available anywhere. give an error
    # or return NULL
    if(required) {
      stop(paste0("need column '", par, "' in data_daily or parameter '", par_, "' in specs"))
    } else {
      return(NULL)
    }
  }
  
  # create the fit as an empty data.frame of the same nrow as data_daily
  fit <- select(data_daily, date)
  
  # add n, the nrow of fit and data_daily
  assign('n', nrow(fit), envir=env)
  
  # add discharge.daily to environment if possible. It's required if pool_K600 
  # uses discharge.daily; otherwise, still try (but with error tolerance) in
  # case GPP_daily, etc. are calculated from it by user function
  assign('discharge.daily', get_par('discharge.daily', required=features$pool_K600 %in% c('linear','binned')), envir=env)
  
  # support hierarchical simulation if requested
  if(features$pool_K600 == 'binned') {
    assign('K600_lnQ_nodes_centers', get_par('K600_lnQ_nodes_centers'), envir=env)
    assign('lnK600_lnQ_nodes', get_par('lnK600_lnQ_nodes'), envir=env)
    assign('lnK600_daily_predlog', predict_Kb(K600_lnQ_nodes_centers, lnK600_lnQ_nodes, log(discharge.daily)), envir=env)
  }  
  
  # get daily parameter needs
  needs <- unlist(get_param_names(metab_model)) # keep names to know whether required
  # sort needs to match data_daily default order, which is the order of
  # operations we want to support
  ops.order <- names(eval(formals('metab_sim')$data_daily))
  needs <- needs[na.omit(match(ops.order, needs))]
  
  # add columns to data_daily for each recognized need
  for(needname in names(needs)) {
    need <- needs[needname]
    fit[need] <- get_par(need, required=grepl('^required', needname))
    assign(need, fit[[need]], envir=env)
  }
  
  # package data_daily as the 'fitted' parameters for this simulation run
  metab_model@fit <- fit
  
  # use default get_params() to package output more nicely
  NextMethod()
}

#' @describeIn predict_DO Simulate values for DO.obs (with process and 
#'   observation error), DO.mod (with process error only), and DO.pure (with no 
#'   error). The errors are randomly generated on every new call to predict_DO.
#' @export
#' @importFrom stats rnorm
predict_DO.metab_sim <- function(metab_model, date_start=NA, date_end=NA, ...) {
  
  specs <- get_specs(metab_model)
  sim_seed <- specs$sim_seed
  if(!is.na(sim_seed)) set.seed(sim_seed)
  
  # simulate the daily parameters and fix them for the remainder of this function call
  metab_model@fit <- get_params(metab_model)
  
  # simulate errors to add to modeled data
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

#' Simulate a lnK ~ lnQ relationship in the Kb format
#' 
#' Uses linear interpolation among "nodes" (lnQ, lnK points) to describe the lnK
#' ~ lnQ relationship
#' 
#' @param specs list, preferably generated by specs(), that includes 
#'   K600_lnQ_cnode_meanlog, K600_lnQ_cnode_sdlog, K600_lnQ_nodediffs_meanlog,
#'   and K600_lnQ_nodediffs_sdlog
#' @param K600_lnQ_nodes_centers numeric vector of centers that replaces the 
#'   same element in specs
#' @export
sim_Kb <- function(specs, K600_lnQ_nodes_centers) {
  
  # revise/add specs$K600_lnQ_nodes_centers if second argument was provided
  if(!missing(K600_lnQ_nodes_centers)) specs$K600_lnQ_nodes_centers <- K600_lnQ_nodes_centers
  
  # check Q bins
  if(!is.numeric(K600_lnQ_nodes_centers))
    stop("K600_lnQ_nodes_centers must be numeric")
  
  # simulate piecewise (binned) relationship for lnK600 ~ lnQ
  lnK600_lnQ_nodes <- numeric(length(K600_lnQ_nodes_centers)) # initialize vector
  cnode <- ceiling(length(K600_lnQ_nodes_centers) / 2) # central node, or node 0.5 past center
  lnK600_lnQ_nodes[cnode] <- rnorm(1, K600_lnQ_cnode_meanlog, K600_lnQ_cnode_sdlog)
  if(cnode > 1)
    for(i in rev(seq_len(cnode-1))) {
      lnK600_lnQ_nodes[i] <- rnorm(1, lnK600_lnQ_nodes[i+1] - K600_lnQ_nodediffs_meanlog, K600_lnQ_nodediffs_sdlog)
    }
  if(length(K600_lnQ_nodes_centers) > cnode) 
    for(i in (cnode+1):length(K600_lnQ_nodes_centers)) {
      lnK600_lnQ_nodes[i] <- rnorm(1, lnK600_lnQ_nodes[i-1] + K600_lnQ_nodediffs_meanlog, K600_lnQ_nodediffs_sdlog)
    }
  lnK600_lnQ_nodes
}

#' Predict ln(K600) as Kb-style function of discharge
#' 
#' Uses linear interpolation among "nodes" (lnQ, lnK points) to predict daily 
#' values of the natural log of K600, based on the lnK ~ lnQ relationship 
#' specified by \code{K600_lnQ_nodes_centers} and \code{lnK600_lnQ_nodes}
#' 
#' @inheritParams specs
#' @param lnQ.daily vector of daily values of the natural log of discharge, 
#'   e.g., \code{log(data_daily$discharge.daily)}
#' @export
predict_Kb <- function(K600_lnQ_nodes_centers, lnK600_lnQ_nodes, lnQ.daily) {
  
  # this function is HIGHLY REDUNDANT with metab_bayes.R. See GH#236
  
  # Convert lnQ.daily to bins and bin weights suitable for linear interpolation 
  # from node to node, horizontal at the edges
  bounds <- c(-Inf, K600_lnQ_nodes_centers, Inf)
  cuts <- cut(lnQ.daily, breaks=bounds, ordered_result=TRUE)
  widths <- diff(bounds)[cuts]
  bins <- rbind(pmax(1, as.numeric(cuts) - 1), pmin(length(K600_lnQ_nodes_centers), as.numeric(cuts)))
  weights <- ifelse(is.infinite(widths), 1, (bounds[as.numeric(cuts)+1] - lnQ.daily)/widths)
  # package info
  lnQ.bin1 = bins[1,]
  lnQ.bin2 = bins[2,]
  lnQ.bin1.weight = weights
  lnQ.bin2.weight = 1-weights
  
  # Predict lnK600_daily_predlog (lnK600 for each value of discharge.daily)
  lnK600_lnQ_nodes[lnQ.bin1] * lnQ.bin1.weight + lnK600_lnQ_nodes[lnQ.bin2] * lnQ.bin2.weight
}

