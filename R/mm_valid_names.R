#' Get the valid names for a given model type or types
#' 
#' Returns a vector of the \code{model_name}s for the type[s] indicated. If 
#' \code{type} is not supplied, all model types will be included. After being 
#' returned from this function, model names may be translated to something 
#' slightly more readable with \code{\link{mm_parse_name}} if desired.
#' 
#' @inheritParams mm_name
#' @examples
#' mm_valid_names('mle')
#' @export
mm_valid_names <- function(type=c('bayes','mle','night','Kmodel','sim')) {

  type <- match.arg(type, several.ok=TRUE)
  
  # if more than one type is supplied or implied, return valid names for all 
  # types
  if(length(type) > 1) {
    return(do.call(c, lapply(type, mm_valid_names)))
  }
  
  # if just one type is supplied, determine the list of acceptable names. method
  # differs by model type
  switch(
    type,
    bayes={
      # get the list of prepared mcmc files from the 'models' directory. this
      # line is why mm_generate_mcmc_file can't call this function (via
      # mm_names(check_validity=TRUE))
      mnames <- grep('^b_', dir(system.file("models", package="streamMetabolizer")), value=TRUE)
      # rename so our favorite is first
      favorites <- c("b_np_oipi_pm_km.stan","b_np_oipi_pm_km.jags")
      c(favorites, mnames[-which(mnames %in% favorites)])
    },
    mle={
      opts <- expand.grid(
        type='mle',
        pool_K600=c('none'),
        err_obs_iid=c(TRUE, FALSE),
        err_proc_acor=FALSE,
        err_proc_iid=c(FALSE, TRUE),
        ode_method=c('pairmeans','Euler'),
        deficit_src='DO_mod',
        engine=c('nlm'),
        stringsAsFactors=FALSE)
      incompatible <- (opts$err_obs_iid == opts$err_proc_iid)
      opts <- opts[!incompatible, ]
      mnames <- sapply(seq_len(nrow(opts)), function(i) do.call(mm_name, c(opts[i,], list(check_validity=FALSE))))
      # rename so our favorite is first
      favorites <- c("m_np_oi_pm_km.nlm")
      c(favorites, mnames[-which(mnames %in% favorites)])
    },
    night=c(
      # this causes finite recursion because all args are specified and check_validity=FALSE, so mm_name doesn't call mm_valid_names
      mm_name(type='night', pool_K600='none', err_obs_iid=FALSE, err_proc_acor=FALSE, err_proc_iid=TRUE, ode_method="Euler", deficit_src='NA', engine='lm', check_validity=FALSE)
    ),
    Kmodel=c(
      # this causes finite recursion because all [Kmodel] args are specified and check_validity=FALSE, so mm_name doesn't call mm_valid_names
      mm_name(type='Kmodel', engine='lm', check_validity=FALSE),
      mm_name(type='Kmodel', engine='mean', check_validity=FALSE),
      mm_name(type='Kmodel', engine='loess', check_validity=FALSE)
    ),
    sim=c(
      # this causes finite recursion because all args are specified and check_validity=FALSE, so mm_name doesn't call mm_valid_names
      mm_name(type='sim', pool_K600='none', err_obs_iid=TRUE, err_proc_acor=TRUE, err_proc_iid=TRUE, ode_method="pairmeans", deficit_src='NA', engine='rnorm', check_validity=FALSE),
      mm_name(type='sim', pool_K600='none', err_obs_iid=TRUE, err_proc_acor=TRUE, err_proc_iid=TRUE, ode_method="Euler", deficit_src='NA', engine='rnorm', check_validity=FALSE)
    )
  )
}

