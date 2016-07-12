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
#' NLL(metab.pars=c(GPP.daily=4, ER.daily=-7, K600.daily=15))
#' NLL(metab.pars=c(GPP.daily=2, ER.daily=-2, K600.daily=25))
#' nlm(NLL, p=c(GPP.daily=2, ER.daily=-2, K600.daily=25))
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
  
  function(metab.pars) {
    # use numerical integration to predict the timeseries of DO.mod
    DO.mod <- calc_DO(setNames(as.list(metab.pars), par.names))
    
    # calculate & return the negative log likelihood of mod values relative to
    # obs values. equivalent to Bob's original code & formula at
    # http://www.statlect.com/normal_distribution_maximum_likelihood.htm
    if(err_obs_iid) {
      diffs <- DO.obs - DO.mod
      sd <- if('err_obs_iid_sigma' %in% metab.pars) metab.pars$err_obs_iid_sigma else sqrt(mean(diffs^2))
    } else if(err_proc_iid) {
      dDOdt.mod <- diff(DO.mod)
      diffs <- dDOdt.obs - dDOdt.mod
      sd <- if('err_proc_iid_sigma' %in% metab.pars) metab.pars$err_proc_iid_sigma else sqrt(mean(diffs^2))
    }
    -sum(dnorm(x=diffs, sd=sd, log=TRUE))
  }
}
