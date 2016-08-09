#' Create a function to compute the negative log likelihood of a set of
#' metabolism parameter values
#'
#' @param calc_dDOdt a function as from \code{create_calc_dDOdt}
#' @inheritParams mm_name
#' @return a function that will return a negative log likelihood of the data
#'   given a set of metab.pars
#' @importFrom stats dnorm
#' @examples
#' \dontrun{
#' data <- data_metab('1','30')[seq(1,48,by=2),]
#' dDOdt <- create_calc_dDOdt(data, ode_method='trapezoid', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, err_obs_iid=TRUE)
#' NLL <- create_calc_NLL(DO)
#' iter=10; NLL(metab.pars=c(GPP.daily=2, ER.daily=-2, K600.daily=25))
#' iter=20; NLL(metab.pars=c(GPP.daily=4, ER.daily=-7, K600.daily=15))
#' NLL2 <- create_calc_NLL(DO, par.names=c('GPP.daily','ER.daily'))
#' iter=30; NLL2(metab.pars=c(GPP.daily=2, ER.daily=-2), K600.daily=25)
#' iter=NA; nlm(NLL, p=c(GPP.daily=2, ER.daily=-2, K600.daily=25))
#' iter=NA; nlm(NLL2, p=c(GPP.daily=2, ER.daily=-2), K600.daily=31.265)
#' }
#' @export
create_calc_NLL <- function(
  calc_DO,
  par.names=environment(environment(calc_DO)$calc_dDOdt)$metab.needs,
  err_obs_iid=environment(calc_DO)$err_obs_iid,
  err_proc_iid=environment(calc_DO)$err_proc_iid) {
  if(!xor(err_obs_iid, err_proc_iid))
    stop("need err_obs_iid or err_proc_iid but not both or neither")
  
  # pull out info from the calc_dDOdt closure
  DO.obs <- environment(calc_DO)$DO.obs
  
  # pre-calculate anything we can
  if(err_proc_iid) {
    dDOdt.obs <- diff(DO.obs)
  }
  
  function(metab.pars, ...) {
    # use numerical integration to predict the timeseries of DO.mod. Assumes 
    all.pars <- c(setNames(as.list(metab.pars), par.names), list(...))
    DO.mod <- calc_DO(all.pars)
    
    # calculate & return the negative log likelihood of mod values relative to
    # obs values. equivalent to Bob's original code & formula at
    # http://www.statlect.com/normal_distribution_maximum_likelihood.htm
    if(err_obs_iid) {
      diffs <- DO.obs - DO.mod
      sd <- if('err_obs_iid_sigma' %in% all.pars) all.pars$err_obs_iid_sigma else sqrt(mean(diffs^2))
    } else if(err_proc_iid) {
      dDOdt.mod <- diff(DO.mod) # should this be diff(DO.mod) or calc_dDOdt over all times?
      diffs <- dDOdt.obs - dDOdt.mod
      sd <- if('err_proc_iid_sigma' %in% all.pars) all.pars$err_proc_iid_sigma else sqrt(mean(diffs^2))
    }
    nll <- -sum(dnorm(x=diffs, sd=sd, log=TRUE))
    it <- get('iter', envir=.GlobalEnv)
    if(is.na(it)) {
      iter <<- it <- 1
      plot(DO.mod, type='l', ylim=c(0,20))
      points(DO.obs, pch=19, col='blue', cex=0.5)
    }
    iter <<- it + 1
    lines(DO.mod, col=grey(0.95-min(0.95,it/150)))
    nll
  }
}
