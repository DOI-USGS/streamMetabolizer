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
#' @param index integer or logical. Applicable only if plotting posteriors, and 
#'   useful only if the parname is for a parameter having multiple (e.g., daily)
#'   instances. In this case, the index selects the instance and corresponds to 
#'   the row number in the data.frame element of \code{get_fit(metab_model)} 
#'   that contains the parameter, e.g. \code{get_fit(metab_model)$daily} for 
#'   \code{'GPP_daily'}. The default, TRUE, selects and pools all instances of
#'   the parameter.
#' @param style character indicating which graphics package to use
#' @import dplyr
#' @importFrom tidyr spread
#' @importFrom stats dunif qnorm dnorm qlnorm dlnorm qbeta dbeta qgamma dgamma
#' @export
#' @examples
#' mm_priors_only <- metab_model(specs=specs('bayes', K600_daily_mu=30))
#' plot_distribs(mm_priors_only, 'err_proc_iid_sigma')
plot_distribs <- function(
  metab_model, 
  parname=c('GPP_daily','ER_daily','K600_daily','K600_daily_mu','K600_daily_beta','K600_daily_sigma',
            'err_obs_iid_sigma','err_proc_acor_phi','err_proc_acor_sigma','err_proc_iid_sigma'), 
  index=TRUE,
  style=c('dygraphs','ggplot2')) {
  
  style <- match.arg(style)
  if(!class(metab_model)[1] %in% c('metab_model','metab_bayes'))
    stop("can only plot distribs for models of class 'metab_model' (only) or 'metab_bayes'")
  
  # extract just the parameters we're interested in. do it this way rather than
  # with match.arg because we can provide a more useful error message here
  sp <- get_specs(metab_model)
  hpspecs <- if(length(parname) != 1) c() else which(grepl(paste0("^", parname, "_"), names(sp)))
  knownpars <- eval(formals(plot_distribs)$parname)
  if(length(hpspecs) == 0 || !parname %in% knownpars) {
    couldabeen <- unlist(lapply(knownpars, function(kp) if(any(grepl(paste0("^", kp, "_"), names(sp)))) kp else c()))
    stop("could not find ", parname, " hyperparameters in get_specs(metab_model). try one of these: ", paste0(couldabeen, collapse=', '))
  }
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
  
  # create a data.frame illustrating the prior distribution
  densdf <- switch(
    distrib,
    uniform={
      xlim_prior <- c(min=hyperpars$min, max=hyperpars$max)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        dist='prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dunif(x, min=hyperpars$min, max=hyperpars$max))
    },
    normal={
      xlim_prior <- qnorm(c(0.001, 0.999), mean=hyperpars$mu, sd=hyperpars$sigma)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        dist='prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dnorm(x, mean=hyperpars$mu, sd=hyperpars$sigma))
    },
    lognormal={
      xlim_prior <- qlnorm(c(0.0001, 0.9), meanlog=hyperpars$location, sdlog=hyperpars$scale)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        dist='prior',
        x = exp(seq(log(xlim[1]), log(xlim[2]), length.out=1000)),
        y = dlnorm(x, meanlog=hyperpars$location, sdlog=hyperpars$scale))
    },
    beta={
      xlim_prior <- qbeta(c(0, 1), shape1=hyperpars$alpha, shape2=hyperpars$beta)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        dist='prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dbeta(x, shape1=hyperpars$alpha, shape2=hyperpars$beta))
    },
    gamma={
      xlim_prior <- qgamma(c(0, 0.99), shape=hyperpars$shape, rate=hyperpars$rate)
      xlim <- xlim_prior # eventually join this with xlim posterior if available
      data_frame(
        dist='prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dgamma(x, shape=hyperpars$shape, rate=hyperpars$rate))
    },
    stop('unrecognized distribution function'))
  
  # if available, augment the data.frame with data illustrating the posterior
  plot_posterior <- class(metab_model)[1] == 'metab_bayes' && !is.null(get_mcmc(metab_model))
  if(plot_posterior) {
    mc <- get_mcmc(metab_model)
    # extract MCMC draws from the specified index/indices. If indices, collapse
    # from matrix into vector
    draws <- extract(mc, pars=parname)[[parname]]
    indexed <- is.matrix(draws)
    if(indexed) draws <- c(draws[,index])
    # generate density w/ 1000 points along the line
    post <- density(draws, n=1000)[c('x','y')] %>% 
      as_data_frame() %>%
      mutate(dist='posterior') %>%
      select(dist, x, y)
    densdf <- bind_rows(densdf, post)
    if(!indexed && !missing(index)) warning('index will be ignored because posterior is not indexed')
  } else {
    indexed <- FALSE
    if(!missing(index)) warning('index will be ignored because priors are never indexed and posterior is unavailable')
  }
  
  # prepare the plot title
  ptitle <- paste0(
    parname, 
    if(indexed) paste0(
      '[',
      if(is.logical(index) || length(index) == 1) { 
        as.character(index)
      } else if(all(diff(index) == 1)) {
        paste0(min(index),':',max(index))
      } else {
        paste0('c(', paste0(index, collapse=','), ')')
      },
      ']'))
  
  # make the plot
  p <- switch(
    style,
    ggplot2={
      if(!requireNamespace("ggplot2", quietly=TRUE))
        stop("call install.packages('ggplot2') before plotting with style='ggplot2'")
      dist <- x <- y <- '.ggplot2Var'
      g <- ggplot2::ggplot(densdf, ggplot2::aes(x=x, y=y, fill=dist, color=dist)) + 
        ggplot2::scale_color_manual(values=c(prior='red', posterior='blue')) +
        ggplot2::scale_fill_manual(values=c(prior='red', posterior='blue')) +
        ggplot2::geom_area(alpha=0.4) + 
        ggplot2::theme_bw() +
        ggplot2::ggtitle(ptitle)
      suppressWarnings(print(g))
      g
    },
    dygraphs={
      if(!requireNamespace("dygraphs", quietly=TRUE))
        stop("call install.packages('dygraphs') before plotting with style='dygraphs'")
      # prepare the data for dygraphs. if the two distributions overlap on the x
      # axis, they'll look really funny unless we fill in the NA values, so also
      # approx those in
      dydensdf <- spread(densdf, dist, y)
      if(plot_posterior) {
        dydensdf <- dydensdf %>%
          arrange(x) %>%
          mutate(
            prior = approx(x=x, y=prior, xout=x)$y,
            posterior = approx(x=x, y=posterior, xout=x)$y
          )
      }
      d <- dygraphs::dygraph(dydensdf, main=ptitle) %>% 
        dygraphs::dyAxis('x', rangePad = 5) %>%
        dygraphs::dyOptions(fillAlpha = 0.4) %>%
        dygraphs::dyRangeSelector(height = 20) %>% 
        dygraphs::dySeries('prior', color='red', fillGraph = TRUE)
      if(plot_posterior) {
        d <- d %>% dygraphs::dySeries('posterior', color='blue', fillGraph = TRUE)
      }
      print(d)
      d
    })
  
  invisible(p)
}
