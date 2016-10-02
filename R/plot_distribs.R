#' Plot the prior/posterior distributions of a parameter
#' 
#' Plot the prior and/or posterior disitrubtions as implied by the 
#' hyperparameters in a specs list and/or the
#' 
#' @param metab_model a metab_model object from which the specs (for priors) and
#'   parameters (for posteriors) will be derived. To plot priors only, insert a 
#'   complete specs list into a shell of a metab_model with 
#'   \code{metab_model(specs=sp)} and assign the result to this argument
#' @param parname character. the name of the parameter whose distribution[s] you
#'   wish to plot
#' @import dplyr
#' @importFrom stats rnorm runif rgamma rlnorm
#' @export
#' @examples
#' mm_priors_only <- metab_model(specs=specs('bayes', K600_daily_mu=30))
#' plot_distribs(mm_priors_only, 'err_proc_iid_sigma')
plot_distribs <- function(
  metab_model, 
  parname=c('GPP_daily','ER_daily','K600_daily','K600_daily_mu','K600_daily_beta','K600_daily_sigma',
            'err_obs_iid_sigma','err_proc_acor_phi','err_proc_acor_sigma','err_proc_iid_sigma'), 
  style=c('dygraphs','ggplot2')) {
  
  parname <- match.arg(parname)
  style <- match.arg(style)
  if(!class(metab_model)[1] %in% c('metab_model','metab_bayes'))
    stop("can only plot distribs for models of class 'metab_model' (only) or 'metab_bayes'")
  
  # extract just the parameters we're interested in
  sp <- get_specs(metab_model)
  hpspecs <- grepl(paste0("^", parname), names(sp))
  if(length(hpspecs) == 0) stop("could not find ", parname, " hyperparameters in get_specs(metab_model)")
  hyperpars <- sp[hpspecs]
  names(hyperpars) <- substring(names(hyperpars), nchar(parname) + 2)
  
  # determine the appropriate distribution
  distrib <- c(
    GPP_daily='normal', 
    ER_daily='normal',
    K600_daily='normal',
    K600_daily_mu='normal',
    K600_daily_beta='normal',
    K600_daily_sigma='lognormal',
    err_obs_iid_sigma='lognormal',
    err_proc_acor_phi='beta',
    err_proc_acor_sigma='lognormal',
    err_proc_iid_sigma='lognormal'
  )[parname]
  
  densdf <- switch(
    distrib,
    uniform={
      xlim_prior <- c(min=hyperpars$min, max=hyperpars$max)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        x = seq(xlim[1], xlim[2], length.out=1000),
        prior = dunif(x, min=hyperpars$min, max=hyperpars$max))
    },
    normal={
      xlim_prior <- qnorm(c(0.001, 0.999), mean=hyperpars$mu, sd=hyperpars$sigma)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        x = seq(xlim[1], xlim[2], length.out=1000),
        prior = dnorm(x, mean=hyperpars$mu, sd=hyperpars$sigma))
    },
    lognormal={
      xlim_prior <- qlnorm(c(0.0001, 0.9), meanlog=hyperpars$location, sdlog=hyperpars$scale)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        x = exp(seq(log(xlim[1]), log(xlim[2]), length.out=1000)),
        prior = dlnorm(x, meanlog=hyperpars$location, sdlog=hyperpars$scale))
    },
    beta={
      xlim_prior <- qbeta(c(0, 1), shape1=hyperpars$alpha, shape2=hyperpars$beta)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        x = seq(xlim[1], xlim[2], length.out=1000),
        prior = dbeta(x, shape1=hyperpars$alpha, shape2=hyperpars$beta))
    },
    gamma={
      xlim_prior <- qgamma(c(0, 0.99), shape=hyperpars$shape, rate=hyperpars$rate)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        x = seq(xlim[1], xlim[2], length.out=1000),
        prior = dgamma(x, shape=hyperpars$shape, rate=hyperpars$rate))
    },
    stop('unrecognized distribution function'))
  
  p <- switch(
    style,
    ggplot2={
      if(!requireNamespace("ggplot2", quietly=TRUE))
        stop("call install.packages('ggplot2') before plotting with style='ggplot2'")
      ggplot(densdf, aes(x=x, y=prior)) + geom_area(fill='blue', color='blue', alpha=0.4) + theme_bw()
    },
    dygraphs={
      if(!requireNamespace("dygraphs", quietly=TRUE))
        stop("call install.packages('dygraphs') before plotting with style='dygraphs'")
      dygraph(densdf) %>% dySeries('prior', color='blue', fillGraph = TRUE) %>% 
        dyOptions(fillAlpha = 0.4) %>%
        dygraphs::dyRangeSelector(height = 20)
    })
  
  print(p)
  invisible(p)
}
