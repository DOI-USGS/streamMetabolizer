#' Create a function to compute the numerical integration of a dDOdt function
#' 
#' @param calc_dDOdt a function as from \code{create_calc_dDOdt}
#' @inheritParams mm_name
#' @param err.obs optional numerical vector of length nrow(data) in units of gO2
#'   m^3. Appropriate for simulation, when this vector of observation errors
#'   will be added to the calculated DO values to simulate observation error.
#'   But usually (for MLE or prediction from a fitted MLE/Bayesian/nighttime
#'   regression model) \code{err.obs} should be missing or 0
#' @return a function that will return a negative log likelihood of the data 
#'   given a set of metab.pars
#' @examples
#' \dontrun{
#' # prepare data for examples
#' data <- data_metab('3','30')[97:144,][seq(1,48,by=2),]
#' # preds.init <- list(GPP.daily=2.82,ER.daily=-2.12,K600.daily=31.27)
#' preds.init <- as.list(dplyr::select(
#'   predict_metab(metab(specs(mm_name('mle', ode_method='trapezoid')), data=data)),
#'   GPP.daily=GPP, ER.daily=ER, K600.daily=K600))
#' DOtime <- data$solar.time
#' dDOtime <- data$solar.time[-nrow(data)] + (data$solar.time[2] - data$solar.time[1])/2
#' 
#' # integration of dDOdt by euler, trapezoid, rk2, rk4, and lsoda methods
#' plot(x=DOtime, y=data$DO.obs, pch=3, cex=0.6)
#' # euler
#' dDOdt <- create_calc_dDOdt(data, ode_method='euler', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='euler')
#' DO.mod.euler <- DO(metab.pars=preds.init)
#' lines(x=DOtime, y=DO.mod.euler, type='l', col='chartreuse3')
#' # trapezoid=pairmeans
#' dDOdt <- create_calc_dDOdt(data, ode_method='trapezoid', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='trapezoid')
#' DO.mod.trap <- DO(metab.pars=preds.init)
#' lines(x=DOtime, y=DO.mod.trap, type='l', col='gold')
#' # lsoda
#' dDOdt <- create_calc_dDOdt(data, ode_method='lsoda', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='lsoda')
#' DO.mod <- DO(metab.pars=preds.init)
#' lines(x=DOtime, y=DO.mod, type='l', col='navy')
#' # rk2
#' dDOdt <- create_calc_dDOdt(data, ode_method='rk2', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='rk2')
#' DO.mod <- DO(metab.pars=preds.init)
#' lines(x=DOtime, y=DO.mod, type='l', col='blue')
#' # rk4
#' dDOdt <- create_calc_dDOdt(data, ode_method='rk4', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='rk4')
#' DO.mod <- DO(metab.pars=preds.init)
#' lines(x=DOtime, y=DO.mod, type='l', col='magenta')
#' 
#' # with observation and/or process error
#' plot(x=DOtime, y=data$DO.obs, col='black', pch=3, cex=0.6)
#' dDOdt <- create_calc_dDOdt(data, ode_method='trapezoid', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='trapezoid', err.obs=rnorm(nrow(data), 0, 0))
#' DO.mod.obs <- DO(preds.init)
#' lines(x=DOtime, y=DO.mod.obs, col='black', lty=1)
#' dDOdt <- create_calc_dDOdt(data, ode_method='trapezoid', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='trapezoid', err.obs=rnorm(nrow(data), 0, 0.3))
#' DO.mod.obserr <- DO(preds.init)
#' lines(x=DOtime, y=DO.mod.obserr, col='blue', lty=2)
#' dDOdt <- create_calc_dDOdt(data, ode_method='trapezoid', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod', err.proc=rnorm(nrow(data), 0, 2))
#' DO <- create_calc_DO(dDOdt, ode_method='trapezoid', err.obs=rnorm(nrow(data), 0, 0.1))
#' DO.mod.operr <- DO(preds.init)
#' lines(x=DOtime, y=DO.mod.operr, col='red', lty=2)
#' }
#' @export
create_calc_DO <- function(calc_dDOdt, ode_method=environment(calc_dDOdt)$ode_method, err.obs=0) {
  
  # pull out info from the calc_dDOdt closure (all from data)
  DO.obs.1 <- environment(calc_dDOdt)$data$DO.obs[1]
  t <- environment(calc_dDOdt)$data$t
  
  if(requireNamespace('deSolve', quietly=TRUE) && !(ode_method %in% c('Euler','pairmeans'))) {
    # identify the right ode method argument
    ode.method <- switch(
      ode_method,
      euler=, trapezoid='euler', # we do the trapezoidy/pairmeansy stuff in calc_dDOdt
      rk2=deSolve::rkMethod('rk2'),
      ode_method
    )
    
    # use numerical integration to predict the timeseries of DO.mod
    calc.DO <- function(metab.pars) {
      DO.mod.1 <- if(exists('DO.mod.1', metab.pars)) metab.pars$DO.mod.1 else DO.obs.1
      DO.mod <- deSolve::ode(
        y=c(DO.mod=DO.mod.1),
        parms=metab.pars,
        times=t,
        func=calc_dDOdt, method=ode.method)[,'DO.mod']
      DO.mod + err.obs
    }
  } else {
    # identify the right ode method argument
    ode.method <- switch(
      ode_method,
      Euler=, pairmeans='euler', # we do the trapezoidy/pairmeansy stuff in calc_dDOdt
      stop("package deSolve is required for ode_method '", ode_method, "'.\n",
           "  Either install deSolve or select ode_method from c('euler','trapezoid')")
    )
    
    # use numerical integration to predict the timeseries of DO.mod
    calc.DO <- function(metab.pars) {
      DO.mod.1 <- if(exists('DO.mod.1', metab.pars)) metab.pars[['DO.mod.1']] else DO.obs.1
      DO.mod <- c(DO.mod.1, rep(NA, length(t)-1))
      for(i in t[-1]) {
        DO.mod[i] <-
          DO.mod[i-1] +
          calc_dDOdt(t=i-1, state=c(DO.mod=DO.mod[i-1]), metab.pars=metab.pars)$dDOdt
      }
      DO.mod + err.obs
    }
  }
  calc.DO
}
