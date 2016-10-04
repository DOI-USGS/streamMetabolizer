#' Plot the prior/posterior distributions of a parameter
#' 
#' Plot the prior and/or posterior disitrubtions as implied by the 
#' hyperparameters in a specs list and/or the
#' 
#' @param dist_data Either a specs list (for priors only) or a metab_model
#'   object (for both priors and posteriors).
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
#' \dontrun
#' # priors only
#' plot_distribs(specs('bayes', K600_daily_mu=30), 'K600_daily')
#' 
#' # posteriors, too
#' mm <- metab(specs(mm_name('bayes')), data=data_metab('1', res='30'))
#' plot_distribs(mm, 'GPP_daily', 1)
#' 
#' # with modifications
#' plot_distribs(mm, 'err_proc_iid_sigma') %>%
#'   dygraphs::dyRangeSelector(dateWindow=c(-0.1,1.3)) %>%
#'   dygraphs::dyAxis(name='y', valueRange=c(0,15))
#' }
plot_distribs <- function(
  dist_data, 
  parname=c('GPP_daily','ER_daily','K600_daily','K600_daily_mu','K600_daily_beta','K600_daily_sigma',
            'err_obs_iid_sigma','err_proc_acor_phi','err_proc_acor_sigma','err_proc_iid_sigma'), 
  index=TRUE,
  style=c('dygraphs','ggplot2'),
  plot_prior_rescaled=TRUE) {
  
  style <- match.arg(style)
  if(!class(dist_data)[1] %in% c('specs','metab_bayes')) {
    stop("can only plot distribs for models of class 'specs' or 'metab_bayes'")
  }
  
  # extract just the parameters we're interested in. do it this way rather than
  # with match.arg because we can provide a more useful error message here
  sp <- switch(
    class(dist_data)[1],
    specs = dist_data,
    metab_bayes = get_specs(dist_data))
  hpspecs <- if(length(parname) != 1) c() else which(grepl(paste0("^", parname, "_"), names(sp)))
  knownpars <- eval(formals(plot_distribs)$parname)
  if(length(hpspecs) == 0 || !parname %in% knownpars) {
    couldabeen <- unlist(lapply(knownpars, function(kp) if(any(grepl(paste0("^", kp, "_"), names(sp)))) kp else c()))
    msg1 <- if(length(parname) > 0) {
      paste0("parname must have length 1.")
    } else {
      paste0("could not find ", parname, " hyperparameters in get_specs(dist_data).")
    }
    stop(msg1, " try one of these: ", paste0(couldabeen, collapse=', '))
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
      xlim <- c(min=hyperpars$min, max=hyperpars$max)
      data_frame(
        dist = 'prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dunif(x, min=hyperpars$min, max=hyperpars$max))
    },
    normal={
      xlim <- qnorm(c(0.001, 0.999), mean=hyperpars$mu, sd=hyperpars$sigma)
      data_frame(
        dist = 'prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dnorm(x, mean=hyperpars$mu, sd=hyperpars$sigma))
    },
    lognormal={
      xlim <- qlnorm(c(0.0001, 0.9), meanlog=hyperpars$location, sdlog=hyperpars$scale)
      prior_rescaled <- exp(hyperpars$location + rnorm(1000000, 0, 1)*hyperpars$scale) %>%
        log() %>% density() %>% # density is smoother if done in log space
        .[c('x','y')] %>%
        as.data.frame() %>%
        mutate(x = exp(x)) # have to stop here and make separate call to mutate(dist=...) to avoid "Error in FUN(left, right) : non-numeric argument to binary operator"
      bind_rows(
        data_frame(
          dist = 'prior',
          x = exp(seq(log(xlim[1]), log(xlim[2]), length.out=1000)),
          y = dlnorm(x, meanlog=hyperpars$location, sdlog=hyperpars$scale)),
        prior_rescaled %>%
          mutate(dist = 'prior_rescaled'))
    },
    beta={
      xlim <- qbeta(c(0, 1), shape1=hyperpars$alpha, shape2=hyperpars$beta)
      data_frame(
        dist = 'prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dbeta(x, shape1=hyperpars$alpha, shape2=hyperpars$beta))
    },
    gamma={
      xlim <- qgamma(c(0, 0.99), shape=hyperpars$shape, rate=hyperpars$rate)
      data_frame(
        dist = 'prior',
        x = seq(xlim[1], xlim[2], length.out=1000),
        y = dgamma(x, shape=hyperpars$shape, rate=hyperpars$rate))
    },
    stop('unrecognized distribution function'))
  plot_prior_rescaled <- plot_prior_rescaled && ('prior_rescaled' %in% densdf$dist)
  if(!plot_prior_rescaled) {
    densdf <- filter(densdf, dist != 'prior_rescaled')
  }
  
  # if available, augment the data.frame with data illustrating the posterior
  plot_posterior <- class(dist_data)[1] == 'metab_bayes' && !is.null(get_mcmc(dist_data))
  if(plot_posterior) {
    mc <- get_mcmc(dist_data)
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
  
  # prepare shared plot style info
  colorvals <- c(prior='red', prior_rescaled='gold', posterior='blue')
  
  # make the plot
  plot_out <- switch(
    style,
    ggplot2={
      if(!requireNamespace("ggplot2", quietly=TRUE))
        stop("call install.packages('ggplot2') before plotting with style='ggplot2'")
      dist <- x <- y <- '.ggplot2Var'
      ggplot2::ggplot(densdf, ggplot2::aes(x=x, y=y, fill=dist, color=dist)) + 
        ggplot2::scale_color_manual(values=colorvals) +
        ggplot2::scale_fill_manual(values=colorvals) +
        ggplot2::geom_area(alpha=0.4, na.rm=TRUE) + 
        ggplot2::theme_bw() +
        ggplot2::ggtitle(ptitle)
    },
    dygraphs={
      if(!requireNamespace("dygraphs", quietly=TRUE))
        stop("call install.packages('dygraphs') before plotting with style='dygraphs'")
      # prepare the data for dygraphs. if the distributions overlap on the x 
      # axis, they'll look really funny unless we fill in the NA values, so also
      # approx those in
      dydensdf <- spread(densdf, dist, y)
      if(plot_prior_rescaled || plot_posterior) {
        dydensdf <- dydensdf %>%
          arrange(x) %>%
          mutate(prior = approx(x=x, y=prior, xout=x)$y)
        if(plot_prior_rescaled)
          dydensdf <- dydensdf %>%
            mutate(prior_rescaled = approx(x=x, y=prior_rescaled, xout=x)$y)
        if(plot_posterior)  
          dydensdf <- dydensdf %>%
            mutate(posterior = approx(x=x, y=posterior, xout=x)$y)
      }
      d <- dygraphs::dygraph(dydensdf, main=ptitle) %>% 
        dygraphs::dyAxis('x', rangePad = 5) %>%
        dygraphs::dyOptions(fillAlpha = 0.4) %>%
        dygraphs::dyRangeSelector(height = 20) %>%
        dygraphs::dySeries('prior', color=colorvals[['prior']], fillGraph = TRUE)
      if(plot_prior_rescaled) {
        d <- d %>% dygraphs::dySeries('prior_rescaled', color=colorvals[['prior_rescaled']], fillGraph = TRUE)
      }
      if(plot_posterior) {
        d <- d %>% dygraphs::dySeries('posterior', color=colorvals[['posterior']], fillGraph = TRUE)
      }
      d
    })
  
  plot_out
}
