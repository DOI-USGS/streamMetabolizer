#' Create a function that generates a 1-day timeseries of DO.mod
#' 
#' Creates a closure that bundles data and helper functions into a single 
#' function that returns dDOdt in gO2 m^-3 timestep^-1 for any given time t.
#' 
#' @section ode_method 'trapezoid':
#' 
#' 'pairmeans' and 'trapezoid' are identical. They are the analytical
#' solution to a trapezoid rule with this starting point:
#' \code{
#'  DO.mod[t+1] = 
#'   DO.mod[t] +
#'    (((GPP[t]+GPP[t+1])/2) / (depth[t]+depth[t+1])/2
#'    + ((ER[t]+ER[t+1])/2) / (depth[t]+depth[t+1])/2 
#'    + (k.O2[t](DO.sat[t] - DO.mod[t]) + k.O2[t+1](DO.sat[t+1] - DO.mod[t+1]))/2
#'    + ((err.proc[t]+err.proc[t+1])/2) / (depth[t]+depth[t+1])/2 
#'    ) * timestep
#' }
#' and this solution:
#' \code{
#' DO.mod[t+1] - DO.mod[t] =
#'  (- DO.mod[t] * (k.O2[t]+k.O2[t+1])/2
#'    + (GPP[t]+GPP[t+1] +
#'       ER[t]+ER[t+1] +
#'       err.proc[t]+err.proc[t+1]
#'      ) / (depth[t]+depth[t+1])
#'    + (k.O2[t]*DO.sat[t] + k.O2[t+1]*DO.sat[t+1])/2 
#'  ) * timestep / (1 + timestep*k.O2[t+1]/2)
#' }
#' where we're treating err.proc as a rate in gO2/m2/d, just like GPP & ER, and
#' err.proc=0 for model fitting.
#' 
#' @param data data.frame as in \code{\link{metab}}, except that data must 
#'   contain exactly one date worth of inputs (~24 hours according to 
#'   \code{\link{specs}$day_start} and \code{\link{specs}$day_end}).
#' @inheritParams mm_name
#' @param err.proc optional numerical vector of length nrow(data). Process error
#'   in units of gO2 m^-2 d^-1 (THIS MAY DIFFER FROM WHAT YOU'RE USED TO!). 
#'   Appropriate for simulation, when this vector of process errors will be 
#'   added to the calculated values of GPP and ER (then divided by depth and
#'   multiplied by timestep duration) to simulate process error. But usually
#'   (for MLE or prediction from a fitted MLE/Bayesian/nighttime regression
#'   model) \code{err.proc} should be missing or 0
#' @return a function that accepts args \code{t} (the time in 0:(n-1) where n is
#'   the number of timesteps), \code{DO.mod.t} (the value of DO.mod at time t in
#'   gO2 m^-3), and \code{metab} (a list of metabolism parameters; to see which 
#'   parameters should be included in this list, create \code{dDOdt} with this 
#'   function and then call \code{environment(dDOdt)$metab.needs})
#' @import dplyr
#' @importFrom unitted v
#' @importFrom stats approxfun
#' @export
#' @examples
#' \dontrun{
#' data <- data_metab('1','30')[seq(1,48,by=2),]
#' dDOdt.obs <- diff(data$DO.obs)
#' preds.init <- as.list(dplyr::select(
#'   predict_metab(metab(specs(mm_name('mle', ode_method='euler')), data=data)),
#'   GPP.daily=GPP, ER.daily=ER, K600.daily=K600))
#' DOtime <- data$solar.time
#' dDOtime <- data$solar.time[-nrow(data)] + (data$solar.time[2] - data$solar.time[1])/2
#' 
#' # args to create_calc_dDOdt determine which values are needed in metab.pars
#' dDOdt <- create_calc_dDOdt(data, ode_method='trapezoid', GPP_fun='satlight',
#'   ER_fun='q10temp', deficit_src='DO_mod')
#' names(formals(dDOdt)) # always the same: args to pass to dDOdt()
#' environment(dDOdt)$metab.needs # get the names to be included in metab.pars
#' dDOdt(t=23, state=c(DO.mod=data$DO.obs[1]),
#'   metab.pars=list(Pmax=0.2, alpha=0.01, ER20=-0.05, K600.daily=3))$dDOdt
#' 
#' # different required args; try in a timeseries
#' dDOdt <- create_calc_dDOdt(data, ode_method='euler', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' environment(dDOdt)$metab.needs # get the names to be included in metab
#' # approximate dDOdt and DO using DO.obs for DO deficits & Eulerian integration
#' DO.mod.m <- data$DO.obs[1]
#' dDOdt.mod.m <- NA
#' for(t in 1:23) {
#'  dDOdt.mod.m[t] <- dDOdt(t=t, state=c(DO.mod=DO.mod.m[t]),
#'     metab.pars=list(GPP.daily=2, ER.daily=-1.4, K600.daily=21))$dDOdt
#'  DO.mod.m[t+1] <- DO.mod.m[t] + dDOdt.mod.m[t]
#' }
#' par(mfrow=c(2,1), mar=c(3,3,1,1)+0.1)
#' plot(x=DOtime, y=data$DO.obs)
#' lines(x=DOtime, y=DO.mod.m, type='l', col='purple')
#' plot(x=dDOtime, y=dDOdt.obs)
#' lines(x=dDOtime, y=dDOdt.mod.m, type='l', col='blue')
#' par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
#' 
#' # compute & plot a full timeseries with ode() integration
#' dDOdt <- create_calc_dDOdt(data, ode_method='euler', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' DO.mod.o <- deSolve::ode(
#'   y=c(DO.mod=data$DO.obs[1]),
#'   parms=list(GPP.daily=2, ER.daily=-1.4, K600.daily=21),
#'   times=1:nrow(data), func=dDOdt, method='euler')[,'DO.mod']
#' par(mfrow=c(2,1), mar=c(3,3,1,1)+0.1)
#' plot(x=DOtime, y=data$DO.obs)
#' lines(x=DOtime, y=DO.mod.m, type='l', col='purple')
#' lines(x=DOtime, y=DO.mod.o, type='l', col='red', lty=2)
#' dDOdt.mod.o <- diff(DO.mod.o)
#' plot(x=dDOtime, y=dDOdt.obs)
#' lines(x=dDOtime, y=dDOdt.mod.m, type='l', col='blue')
#' lines(x=dDOtime, y=dDOdt.mod.o, type='l', col='green', lty=2)
#' par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
#' 
#' # see how values of metab.pars affect the dDOdt predictions
#' library(dplyr); library(ggplot2); library(tidyr)
#' dDOdt <- create_calc_dDOdt(data, ode_method='euler', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod')
#' apply_dDOdt <- function(GPP.daily, ER.daily, K600.daily) {
#'   DO.mod.m <- data$DO.obs[1]
#'   dDOdt.mod.m <- NA
#'   for(t in 1:23) {
#'    dDOdt.mod.m[t] <- dDOdt(t=t, state=c(DO.mod=DO.mod.m[t]),
#'     list(GPP.daily=GPP.daily, ER.daily=ER.daily, K600.daily=K600.daily))$dDOdt
#'    DO.mod.m[t+1] <- DO.mod.m[t] + dDOdt.mod.m[t]
#'   }
#'   dDOdt.mod.m
#' }
#' dDO.preds <- data_frame(
#'   solar.time = dDOtime,
#'   dDO.preds.base = apply_dDOdt(3, -5, 15),
#'   dDO.preds.dblGPP = apply_dDOdt(6, -5, 15),
#'   dDO.preds.dblER = apply_dDOdt(3, -10, 15),
#'   dDO.preds.dblK = apply_dDOdt(3, -5, 30))
#' dDO.preds %>%
#'   gather(key=dDO.series, value=dDO.dt, starts_with('dDO.preds')) %>%
#'   ggplot(aes(x=solar.time, y=dDO.dt, color=dDO.series)) + geom_line() + theme_bw()
#' 
#' # try simulating process eror
#' data <- data_metab('1','30')[seq(1,48,by=2),]
#' plot(x=data$solar.time, y=data$DO.obs)
#' dDOdt.noerr <- create_calc_dDOdt(data, ode_method='rk4', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod', err.proc=rep(0, nrow(data)))
#' DO.mod.noerr <- deSolve::ode(
#'   y=c(DO.mod=data$DO.obs[1]),
#'   parms=list(GPP.daily=2, ER.daily=-1.4, K600.daily=21),
#'   times=1:nrow(data), func=dDOdt.noerr, method='rk4')[,'DO.mod']
#' lines(x=data$solar.time, y=DO.mod.noerr, type='l', col='purple')
#' # with error
#' dDOdt.err <- create_calc_dDOdt(data, ode_method='rk4', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod', err.proc=rep(+0.4, nrow(data)))
#' DO.mod.err <- deSolve::ode(
#'   y=c(DO.mod=data$DO.obs[1]),
#'   parms=list(GPP.daily=2, ER.daily=-1.4, K600.daily=21),
#'   times=1:nrow(data), func=dDOdt.err, method='rk4')[,'DO.mod']
#' lines(x=data$solar.time, y=DO.mod.err, type='l', col='red', lty=2)
#' # with same error each timestep is same as with reduced ER
#' dDOdt.noerr2 <- create_calc_dDOdt(data, ode_method='rk4', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod', err.proc=rep(0, nrow(data)))
#' DO.mod.noerr2 <- deSolve::ode(
#'   y=c(DO.mod=data$DO.obs[1]),
#'   parms=list(GPP.daily=2, ER.daily=-1.4+0.4, K600.daily=21),
#'   times=1:nrow(data), func=dDOdt.noerr2, method='rk4')[,'DO.mod']
#' lines(x=data$solar.time, y=DO.mod.noerr2, type='l', col='green', lty=3)
#' # with different timestep, same error value should mean very similar curve
#' data <- data_metab('1','30')
#' dDOdt.err2 <- create_calc_dDOdt(data, ode_method='rk4', GPP_fun='linlight',
#'   ER_fun='constant', deficit_src='DO_mod', err.proc=rep(0.4, nrow(data)))
#' DO.mod.err2 <- deSolve::ode(
#'   y=c(DO.mod=data$DO.obs[1]),
#'   parms=list(GPP.daily=2, ER.daily=-1.4, K600.daily=21),
#'   times=1:nrow(data), func=dDOdt.err2, method='rk4')[,'DO.mod']
#' lines(x=data$solar.time, y=DO.mod.err2, type='l', col='black', lty=2)
#' }
create_calc_dDOdt <- function(data, ode_method, GPP_fun, ER_fun, deficit_src, err.proc=0) {

  # simplify time indexing. we've guaranteed in mm_model_by_ply that the
  # timesteps are regular
  data$t <- seq_len(nrow(data))

  # define the forcing (temp.water, light, DO.sat, etc.) interpolations and 
  # other inputs to include in the dDOdt() closure. 
  integer.t <- isTRUE(ode_method %in% c('euler','trapezoid','Euler','pairmeans'))
  data$KO2.conv <- convert_k600_to_kGAS(k600=1, temperature=data$temp.water, gas='O2')
  data$err.proc <- err.proc # get replication if needed
  if(integer.t) {
    # for indexing, converting df columns to vectors speeds things up by 40x
    DO.obs <- data$DO.obs
    DO.sat <- data$DO.sat
    temp.water <- data$temp.water
    depth <- data$depth
    light <- data$light
    KO2.conv <- data$KO2.conv
    err.proc <- data$err.proc
  } else { 
    # other methods require functions that can be applied at non-integer 
    # values of t. approxfun is pretty darn fast and ever-so-slightly faster
    # with data$x than with independent vectors of t and a variable
    DO.obs <- approxfun(data$t, data$DO.obs, rule=2)
    DO.sat <- approxfun(data$t, data$DO.sat, rule=2)
    depth <- approxfun(data$t, data$depth, rule=2)
    temp.water <- approxfun(data$t, data$temp.water, rule=2)
    light <- approxfun(data$t, data$light, rule=2)
    KO2.conv <- approxfun(data$t, data$KO2.conv, rule=2)
    err.proc <- if(all(err.proc == 0)) function(t) 0 else approxfun(data$t, data$err.proc, rule=2)
  }
  timestep.days <- suppressWarnings(mean(as.numeric(diff(unitted::v(data$solar.time)), units="days"), na.rm=TRUE))

  # collect the required metab.pars parameter names in a vector called metab.needs
  metab.needs <- c()

  # GPP: instantaneous gross primary production at time t in gO2 m^-2 d^-1
  GPP <- switch(
    GPP_fun,
    'NA'=(function(){
      function(t, metab.pars) 0
    })(), 
    linlight=(function(){
      # normalize light by the sum of light in the first 24 hours of the time window
      mean.light <- with(
        list(in.solar.day = data$solar.time < (data$solar.time[1] + as.difftime(1, units='days'))),
        mean(data$light[in.solar.day]))
      if(mean.light == 0) mean.light <- 1
      metab.needs <<- c(metab.needs, 'GPP.daily')
      if(integer.t) function(t, metab.pars) { 
        metab.pars[['GPP.daily']] * light[t] / mean.light
      } else function(t, metab.pars) {
        metab.pars[['GPP.daily']] * light(t) / mean.light
      }
    })(),
    satlight=(function(){
      metab.needs <<- c(metab.needs, c('Pmax','alpha'))
      if(integer.t) function(t, metab.pars) {
        Pmax <- metab.pars[['Pmax']]; Pmax * tanh(metab.pars[['alpha']] * light[t] / Pmax)
      } else function(t, metab.pars) { 
        Pmax <- metab.pars[['Pmax']]; Pmax * tanh(metab.pars[['alpha']] * light(t) / Pmax)
      }
    })(),
    satlightq10temp=(function(){
      metab.needs <<- c(metab.needs, c('Pmax','alpha'))
      if(integer.t) function(t, metab.pars) {
        Pmax <- metab.pars[['Pmax']]; Pmax * tanh(metab.pars[['alpha']] * light[t] / Pmax) * 1.036 ^ (temp.water[t] - 20)
      } else function(t, metab.pars) {
        Pmax <- metab.pars[['Pmax']]; Pmax * tanh(metab.pars[['alpha']] * light(t) / Pmax) * 1.036 ^ (temp.water(t) - 20)
      }
    })(),
    stop('unrecognized GPP_fun')
  )

  # ER: instantaneous ecosystem respiration at time t in d^-1
  ER <- switch(
    ER_fun,
    constant=(function(){
      metab.needs <<- c(metab.needs, 'ER.daily')
      function(t, metab.pars) {
        metab.pars[['ER.daily']]
      }
    })(),
    q10temp=(function(){
      # song_methods_2016 cite Gulliver & Stefan 1984; Parkhill & Gulliver 1999
      metab.needs <<- c(metab.needs, 'ER20')
      if(integer.t) function(t, metab.pars) {
        metab.pars[['ER20']] * 1.045 ^ (temp.water[t] - 20)
      } else function(t, metab.pars) {
        metab.pars[['ER20']] * 1.045 ^ (temp.water(t) - 20)
      }
    })(),
    stop('unrecognized ER_fun')
  )

  # D: instantaneous reaeration rate at time t in gO2 m^-3 d^-1
  D <- switch(
    deficit_src,
    DO_obs=(function(){
      metab.needs <<- c(metab.needs, 'K600.daily')
      if(integer.t) function(t, metab.pars, DO.mod.t) {
        metab.pars[['K600.daily']] * KO2.conv[t] * (DO.sat[t] - DO.obs[t])
      } else function(t, metab.pars, DO.mod.t) {
        metab.pars[['K600.daily']] * KO2.conv(t) * (DO.sat(t) - DO.obs(t))
      }
    })(),
    DO_obs_filter=, 
    DO_mod=(function(){
      metab.needs <<- c(metab.needs, 'K600.daily')
      if(integer.t) function(t, metab.pars, DO.mod.t) {
        metab.pars[['K600.daily']] * KO2.conv[t] * (DO.sat[t] - DO.mod.t)
      } else function(t, metab.pars, DO.mod.t) {
        metab.pars[['K600.daily']] * KO2.conv(t) * (DO.sat(t) - DO.mod.t)
      }
    })(),
    stop('unrecognized deficit_src')
  )
  
  # dDOdt: instantaneous rate of change in DO at time t in gO2 m^-3 timestep^-1
  dDOdt <- switch(
    ode_method,
    # 'pairmeans' and 'trapezoid' are identical and are the analytical solution
    # to a trapezoid rule. remember we're treating err.proc as a rate in
    # gO2/m2/d, just like GPP & ER
    trapezoid=, pairmeans={
      function(t, state, metab.pars){
        K600.daily <- metab.pars[['K600.daily']]
        list(
          dDOdt={
            - state[['DO.mod']] * K600.daily * {KO2.conv[t]+KO2.conv[t+1]}/2 +
              0.5*{GPP(t, metab.pars) + ER(t, metab.pars) + err.proc[t]}/depth[t] +
              0.5*{GPP(t+1, metab.pars) + ER(t+1, metab.pars) + err.proc[t+1]}/depth[t+1] +
              K600.daily * {KO2.conv[t]*DO.sat[t] + KO2.conv[t+1]*DO.sat[t+1]}/2
          } * timestep.days / {1 + timestep.days * K600.daily * KO2.conv[t+1]/2})
      }
    },
    # all other methods use a straightforward calculation of dDOdt at values of
    # t and DO.mod.t as requested by the ODE solver
    if(integer.t) function(t, state, metab.pars) { # Euler and euler, b/c trapezoid and pairmeans are covered above
      list(
        dDOdt={
          {GPP(t, metab.pars) + ER(t, metab.pars) + err.proc[t]} / depth[t] +
            D(t, metab.pars, state[['DO.mod']])} *
          timestep.days)
    } else function(t, state, metab.pars) {
      list(
        dDOdt={
          {GPP(t, metab.pars) + ER(t, metab.pars) + err.proc(t)} / depth(t) +
            D(t, metab.pars, state[['DO.mod']])} *
          timestep.days)
    }
  )

  # return the closure, which wraps up the final dDOdt code, light(), DO.sat(),
  # depth(), GPP(), ER(), D(), and anything else defined within create_calc_dDOdt
  # into one bundle
  dDOdt
}
