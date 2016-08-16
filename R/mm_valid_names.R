#' Get the valid names for a given model type or types
#' 
#' Returns a vector of the \code{model_name}s for the type[s] indicated. If 
#' \code{type} is not supplied, all model types will be included. After being 
#' returned from this function, model names may be translated to something 
#' slightly more readable with \code{\link{mm_parse_name}} if desired.
#' 
#' @inheritParams mm_name
#' @import dplyr
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
  
  # get lists of all common possibilities
  . <- '.dplyr.var'
  all_ode_methods <- formals(mm_name)$ode_method %>% eval() %>% .[.!='NA']
  all_GPP_funs <- formals(mm_name)$GPP_fun %>% eval() %>% .[.!='NA']
  all_ER_funs <- formals(mm_name)$ER_fun %>% eval() %>% .[.!='NA']
  all_deficit_srcs <- c('DO_mod','DO_obs')
  
  # if just one type is supplied, determine the list of acceptable names. method
  # differs by model type
  mnames <- NA
  switch(
    type,
    bayes={
      # get the list of prepared mcmc files from the 'models' directory. this
      # line is why mm_generate_mcmc_file can't call this function (via
      # mm_names(check_validity=TRUE))
      mnames <- grep('^b_', dir(system.file('models', package='streamMetabolizer')), value=TRUE)
      favorites <- c('b_np_oipi_tr_plrckm.stan','b_np_oipi_tr_plrckm.jags')
    },
    mle={
      opts <- expand.grid(
        type='mle',
        pool_K600='none',
        err_obs_iid=c(TRUE, FALSE),
        err_proc_acor=FALSE,
        err_proc_iid=c(FALSE, TRUE),
        ode_method=all_ode_methods,
        GPP_fun=all_GPP_funs,
        ER_fun=all_ER_funs,
        deficit_src=all_deficit_srcs,
        engine='nlm',
        stringsAsFactors=FALSE)
      incompatible <- (opts$err_obs_iid == opts$err_proc_iid)
      opts <- opts[!incompatible, ]
      favorites <- c("m_np_oi_tr_plrckm.nlm")
    },
    night={
      opts <- expand.grid(
        type='night', 
        pool_K600='none', 
        err_obs_iid=FALSE, 
        err_proc_acor=FALSE, 
        err_proc_iid=TRUE, 
        ode_method='Euler', 
        GPP_fun='NA', 
        ER_fun='constant', 
        deficit_src='DO_obs_filter', 
        engine='lm', 
        stringsAsFactors=FALSE)
      favorites <- c('n_np_pi_eu_rckf.lm')
    },
    Kmodel={
      opts <- expand.grid(
        type='Kmodel',
        engine=c('lm','mean','loess'),
        stringsAsFactors=FALSE)
      favorites <- c('K_Kc___.lm')
    },
    sim={
      opts <- expand.grid(
        type='sim',
        pool_K600='none',
        err_obs_iid=TRUE,
        err_proc_acor=TRUE,
        err_proc_iid=TRUE,
        ode_method=all_ode_methods,
        GPP_fun=all_GPP_funs,
        ER_fun=all_ER_funs,
        deficit_src=all_deficit_srcs,
        engine='rnorm',
        stringsAsFactors=FALSE)
      favorites <- c('s_np_oipcpi_tr_plrckm.rnorm')
      
    }
  )
  
  # create names list if not already done. requires finite recursion because all
  # args are specified and check_validity=FALSE, so mm_name doesn't call
  # mm_valid_names
  if(all(is.na(mnames))) mnames <- sapply(seq_len(nrow(opts)), function(i) suppressWarnings(do.call(mm_name, c(opts[i,], list(check_validity=FALSE)))))
  
  # rename so our favorite is first
  c(favorites, mnames[-which(mnames %in% favorites)])
}

