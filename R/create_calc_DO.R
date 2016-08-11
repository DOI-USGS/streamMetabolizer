#' Create a function to compute the numerical integration of a dDOdt function
#'
#' @param calc_dDOdt a function as from \code{create_calc_dDOdt}
#' @inheritParams mm_name
#' @return a function that will return a negative log likelihood of the data
#'   given a set of metab.pars
#' @examples
#' \dontrun{
#' # prepare data for examples
#' data <- data_metab('3','30')[97:144,][seq(1,48,by=2),]
#' preds.init <- as.list(dplyr::select(
#'   predict_metab(metab(specs(mm_name('mle', ode_method='pairmeans')), data=data)),
#'   GPP.daily=GPP, ER.daily=ER, K600.daily=K600))
#' DOtime <- data$solar.time
#' dDOtime <- data$solar.time[-nrow(data)] + (data$solar.time[2] - data$solar.time[1])/2
#'
#' # integration of dDOdt by Euler, trapezoid, rk2, rk4, and lsoda methods
#' plot(x=DOtime, y=data$DO.obs, pch=3, cex=0.6)
#' # euler
#' dDOdt <- create_calc_dDOdt(data, ode_method='Euler', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='Euler')
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
#' # show that method='euler' really is Euler by several integration implementations
#' plot(x=DOtime, y=data$DO.obs, col='black', pch=3, cex=0.6)
#' # dDOdt
#' dDOdt <- create_calc_dDOdt(data, ode_method='Euler', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO.mod.dDOdt <- data$DO.obs[1]
#' for(t in 2:nrow(data)) { DO.mod.dDOdt[t] <-
#'   DO.mod.dDOdt[t-1] + dDOdt(t-1, c(DO.mod=DO.mod.dDOdt[t-1]), preds.init)$dDOdt }
#' lines(x=DOtime, y=DO.mod.dDOdt, col='purple')
#' # ode
#' DO.mod.ode <- deSolve::ode(y=c(DO.mod=data$DO.obs[1]), parms=preds.init,
#'   times=1:nrow(data), func=dDOdt, method='euler')[,'DO.mod']
#' lines(x=DOtime, y=DO.mod.ode, col='blue')
#' # DO
#' DO <- create_calc_DO(dDOdt, ode_method='Euler')
#' DO.mod.DO <- DO(preds.init)
#' lines(x=DOtime, y=DO.mod.DO, col='chartreuse3')
#' # original calc_DO_mod function
#' DO.mod.old <- do.call(calc_DO_mod, 
#'   c(preds.init, as.list(data[c('DO.sat','depth','temp.water')]), 
#'     list(frac.GPP = data$light/sum(data$light), frac.ER=1/24, frac.D=1/24,
#'     DO.mod.1=data$DO.obs[1], n=nrow(data), ODE_method='Euler')))
#' lines(x=DOtime, y=DO.mod.old, col='black', lty=5)
#'
#' # show that method='trapezoid' really is pairmeans by several implementations
#' plot(x=DOtime, y=data$DO.obs, col='black', pch=3, cex=0.6)
#' # dDOdt
#' dDOdt <- create_calc_dDOdt(data, ode_method='trapezoid', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO.mod.dDOdt <- data$DO.obs[1]
#' for(t in 2:nrow(data)) { DO.mod.dDOdt[t] <-
#'   DO.mod.dDOdt[t-1] + dDOdt(t-1, c(DO.mod=DO.mod.dDOdt[t-1]), preds.init)$dDOdt }
#' lines(x=DOtime, y=DO.mod.dDOdt, col='purple')
#' # DO
#' DO <- create_calc_DO(dDOdt, ode_method='trapezoid')
#' DO.mod.pm.DO <- DO(preds.init)
#' lines(x=DOtime, y=DO.mod.pm.DO, col='forestgreen', lty=3)
#' # rk2 should be similar, but not identical, to pairmeans
#' dDOdt <- create_calc_dDOdt(data, ode_method='rk2', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO <- create_calc_DO(dDOdt, ode_method='rk2')
#' DO.mod.pm.rk2 <- DO(preds.init)
#' lines(x=DOtime, y=DO.mod.pm.rk2, col='red', lty=4)
#' # original calc_DO_mod function
#' DO.mod.old <- do.call(calc_DO_mod, 
#'   c(preds.init, as.list(data[c('DO.sat','depth','temp.water')]), 
#'     list(frac.GPP = data$light/sum(data$light), frac.ER=1/24, frac.D=1/24,
#'     DO.mod.1=data$DO.obs[1], n=nrow(data), ODE_method='pairmeans')))
#' lines(x=DOtime, y=DO.mod.old, col='green', lty=5)
#' }
#' @export
create_calc_DO <- function(calc_dDOdt, ode_method=environment(calc_dDOdt)$ode_method) {
  
  # pull out info from the calc_dDOdt closure
  DO.obs <- environment(calc_dDOdt)$data$DO.obs
  t <- environment(calc_dDOdt)$data$t
  
  if(requireNamespace('deSolve', quietly=TRUE)) {
    # identify the right ode method argument
    ode.method <- switch(
      ode_method,
      Euler=, trapezoid=, pairmeans='euler', # we do the trapezoidy/pairmeansy stuff in calc_dDOdt
      rk2=deSolve::rkMethod('rk2'),
      ode_method
    )
    
    # use numerical integration to predict the timeseries of DO.mod
    calc.DO <- function(metab.pars) {
      DO.mod.1 <- if(exists('DO.mod.1', metab.pars)) metab.pars$DO.mod.1 else environment(calc_dDOdt)$data$DO.obs[1]
      deSolve::ode(
        y=c(DO.mod=DO.mod.1),
        parms=metab.pars,
        times=t,
        func=calc_dDOdt, method=ode.method)[,'DO.mod']
    }
  } else {
    # identify the right ode method argument
    ode.method <- switch(
      ode_method,
      Euler=, trapezoid=, pairmeans='euler', # we do the trapezoidy/pairmeansy stuff in calc_dDOdt
      stop("package deSolve is required for ode_method '", ode_method, "'.\n",
           "  Either install deSolve or select ode_method from c('Euler','trapezoid','pairmeans')")
    )
    
    # use numerical integration to predict the timeseries of DO.mod
    calc.DO <- function(metab.pars) {
      DO.mod.1 <- if(exists('DO.mod.1', metab.pars)) metab.pars$DO.mod.1 else environment(calc_dDOdt)$data$DO.obs[1]
      DO.mod <- c(DO.mod.1, rep(NA, length(t)-1))
      for(i in t[-1]) {
        DO.mod[i] <-
          DO.mod[i-1] +
          calc_dDOdt(t=i-1, state=c(DO.mod=DO.mod[i-1]), metab.pars=metab.pars)$dDOdt
      }
      DO.mod
    }
  }
  calc.DO
}
