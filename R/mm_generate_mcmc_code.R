#' Generate the models in inst/models/bayes
#' 
#' @param bayes_software Which software are we generating code for?
#' @param ode_method The method to use in solving the ordinary differential 
#'   equation for DO. Euler: dDOdt from t=1 to t=2 is solely a function of GPP, 
#'   ER, DO, etc. at t=1. pairmeans: dDOdt from t=1 to t=2 is a function of the 
#'   mean values of GPP, ER, etc. across t=1 and t=2.
#' @param deficit_src From what DO estimate (observed or modeled) should the DO 
#'   deficit be computed?
#' @param err_obs_iid logical. Should IID observation error be included? If not,
#'   the model will be fit to the differences in successive DO measurements,
#'   rather than to the DO measurements themselves.
#' @param err_proc_acor logical. Should autocorrelated process error (with the 
#'   autocorrelation term phi fitted) be included?
#' @param pooling Should the model pool information among days to get more 
#'   consistent daily estimates?
#' @keywords internal
mm_generate_mcmc_file <- function(
  bayes_software=c('jags','stan'), 
  ode_method=c('Euler','pairmeans'),
  deficit_src=c('DO_mod','DO_obs'),
  err_obs_iid=c(TRUE, FALSE),
  err_proc_acor=c(TRUE, FALSE),
  err_proc_iid=c(TRUE, FALSE),
  pooling='none') {
  
  # choose/check arguments
  bayes_software <- match.arg(bayes_software)
  ode_method <- match.arg(ode_method)
  deficit_src <- match.arg(deficit_src)
  err_obs_iid <- if(!is.logical(err_obs_iid)) stop("need err_obs_iid to be a logical of length 1") else err_obs_iid[1]
  err_proc_acor <- if(!is.logical(err_proc_acor)) stop("need err_proc_acor to be a logical of length 1") else err_proc_acor[1]
  err_proc_iid <- if(!is.logical(err_proc_iid)) stop("need err_proc_iid to be a logical of length 1") else err_proc_iid[1]
  pooling <- match.arg(pooling)

  # check argument compatibility
  if(!err_obs_iid && deficit_src == 'DO_mod') stop("if there's no err_obs, deficit_src must be DO_obs")
    
  # name the model
  model_name <- paste0(
    c(none='np', partial='pp')[[pooling]], '_',
    if(err_obs_iid) 'oi', if(err_proc_acor) 'pc', if(err_proc_iid) 'pi', '_',
    c(Euler='eu', pairmeans='pm')[[ode_method]], '_',
    c(DO_mod='km', DO_obs='ko')[[deficit_src]], '.',
    bayes_software
  )
  
  # helper functions
  comment <- function(...) { # prefix with the appropriate comment character[s]
    chr <- switch(bayes_software, jags='#', stan='//')
    paste0(chr, ' ', paste0(...))
  }
  chunk <- function(..., indent=2, newline=TRUE) { # indent a chunk & add a newline
    lines <- c(list(...))
    lines <- lines[!sapply(lines, is.null)]
    lines <- unlist(lines)
    if(newline) lines <- c(lines, '')
    sapply(lines, function(line) {
      paste0(
        paste0(rep(' ',indent),collapse=''), 
        line)
    }, USE.NAMES = FALSE)
  }
  indent <- function(..., indent=2) {
    chunk(..., indent=indent, newline=FALSE)
  }
  f <- function(distrib, ...) { # switch things to JAGS format if needed
    args <- c(list(...))
    if(bayes_software=='jags') {
      switch(
        distrib,
        normal = {
          distrib <- 'dnorm'
          args[2] <- sprintf('pow(%s, -2)', args[2])
        }, 
        uniform = {
          distrib <- 'dunif'
        }
      )
    }
    paste0(
      distrib, '(',
      paste0(args, collapse=', '),
      ')')
  }
  p <- paste0 # partial line: just combine strings into string
  s <- function(...) { # line with stop/semicolon: combine strings into string & add semicolon if needed
    p(p(...), switch(bayes_software, jags='', stan=';'))
  }
  
  # define the model
  model_text <- c(
    
    comment(model_name), '',
    
    ## data ##
    if(bayes_software == 'stan') c(
      
      'data {',
      
      chunk(
        comment('Metabolism distributions'),
        'real GPP_daily_mu;',
        'real GPP_daily_sigma;',
        'real ER_daily_mu;',
        'real ER_daily_sigma;',
        'real K600_daily_mu;',
        'real K600_daily_sigma;'),
      
      chunk(
        comment('Error distributions'),
        if(err_obs_iid) c(
          'real err_obs_iid_sigma_min;',
          'real err_obs_iid_sigma_max;'),
        if(err_proc_acor) c(
          'real err_proc_acor_phi_min;',
          'real err_proc_acor_phi_max;',
          'real err_proc_acor_sigma_min;',
          'real err_proc_acor_sigma_max;'),
        if(err_proc_iid) c(
          'real err_proc_iid_sigma_min;',
          'real err_proc_iid_sigma_max;')),
      
      chunk(
        comment('Daily data'),
        'int <lower=0> n;',
        'real DO_obs_1;'),
      
      indent(
        comment('Data'),
        'vector [n] DO_obs;',
        'vector [n] DO_sat;',
        'vector [n] frac_GPP;',
        'vector [n] frac_ER;',
        'vector [n] frac_D;',
        'vector [n] depth;',
        'vector [n] KO2_conv;'),
      
      '}',''
    ),
    
    ## Stan: transformed data ## - statements evaluated exactly once
    ## JAGS: data ##
    if(bayes_software == 'stan') c(
      'transformed data {'
    ) else if(bayes_software == 'jags') c(
      'data {'
    ),
    
    if(bayes_software == 'stan') c(
      chunk(
        'vector [n-1] coef_GPP;',
        'vector [n-1] coef_ER;',
        switch(
          deficit_src,
          DO_mod = 'vector [n-1] coef_K600_part;',
          DO_obs = 'vector [n-1] coef_K600_full;'
        ),
        if(ode_method == 'pairmeans' && deficit_src == 'DO_mod') 'vector [n-1] DO_sat_pairmean;',
        if(!err_obs_iid) 'vector [n-1] dDO_obs;')
    ),
    
    indent(
      p('for(i in 1:(n-1)) {'),
      indent(
        # Coefficient pre-calculations
        if(ode_method == 'Euler') c(
          comment('Coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)'),
          s('coef_GPP[i]  <- frac_GPP[i] / depth[i]'),
          s('coef_ER[i]   <- frac_ER[ i] / depth[i]'),
          switch(
            deficit_src,
            DO_mod = c(
              s('coef_K600_part[i] <- KO2_conv[i] * frac_D[i]')
            ),
            DO_obs = c(
              p('coef_K600_full[i] <- KO2_conv[i] * frac_D[i] * '),
              s('  (DO_sat[i] - DO_obs[i])')
            )
          )
        ) else if(ode_method == 'pairmeans') c(
          comment('Coefficients by pairmeans (e.g., mean(frac_GPP[i:(i+1)]) applies to the DO step from i to i+1)'),
          s('coef_GPP[i]  <- (frac_GPP[i] + frac_GPP[i+1])/2 / ((depth[i] + depth[i+1])/2)'),
          s('coef_ER[i]   <- (frac_ER[ i] + frac_ER[ i+1])/2 / ((depth[i] + depth[i+1])/2)'),
          switch(
            deficit_src,
            DO_mod = c(
              s('coef_K600_part[i] <- (KO2_conv[i] + KO2_conv[i+1])/2 * (frac_D[i] + frac_D[i+1])/2'),
              s('DO_sat_pairmean[i] <- (DO_sat[i] + DO_sat[i+1])/2')
            ),
            DO_obs = c(
              p('coef_K600_full[i] <- (KO2_conv[i] + KO2_conv[i+1])/2 * (frac_D[i] + frac_D[i+1])/2 *'),
              s('  (DO_sat[i] + DO_sat[i+1] - DO_obs[i] - DO_obs[i+1])/2')
            )
          )
        ),
        
        # dDO pre-calculations
        if(!err_obs_iid) c(
          comment('dDO observations'),
          s('dDO_obs[i] <- DO_obs[i+1] - DO_obs[i]')
        )
      ),
      p('}')
    ),
    
    # close out transformed data (stan) or data (jags)
    c(
      '}',''
    ),
    
    ## Stan: parameters ##
    if(bayes_software == 'stan') c(
      'parameters {',
      
      indent(
        'real GPP_daily;',
        'real ER_daily;',
        'real K600_daily;'),
      
      if(err_proc_acor) indent(
        '',
        'vector [n-1] err_proc_acor_inc;'),
      
      if(err_obs_iid || err_proc_acor || err_proc_iid) indent(
        '',
        if(err_obs_iid) c(
          'real <lower=err_obs_iid_sigma_min,   upper=err_obs_iid_sigma_max>  err_obs_iid_sigma;'),
        if(err_proc_acor) c(
          'real <lower=err_proc_acor_phi_min,   upper=err_proc_acor_phi_max>   err_proc_acor_phi;',
          'real <lower=err_proc_acor_sigma_min, upper=err_proc_acor_sigma_max> err_proc_acor_sigma;'),
        if(err_proc_iid) c(
          'real <lower=err_proc_iid_sigma_min,  upper=err_proc_iid_sigma_max>  err_proc_iid_sigma;')),
      
      '}',''
    ),
    
    ## Stan: transformed parameters ## - statements evaluated once per leapfrog step
    ## JAGS: model ##
    if(bayes_software == 'stan') c(
      'transformed parameters {'
    ) else if(bayes_software == 'jags') c(
      'model {'
    ),
    
    if(bayes_software == 'stan') chunk(
      if(err_obs_iid)
        'vector [n] DO_mod;',
      if(deficit_src == 'DO_obs')
        'vector [n-1] dDO_mod;',
      if(err_proc_acor) 
        'vector [n-1] err_proc_acor;'),

    indent(
      comment('Model DO time series'),
      comment('* ', ode_method,' version'),
      comment('* ', if(!err_obs_iid) 'no ', 'observation error'),
      comment('* ', paste0(c(if(err_proc_iid) 'IID', if(err_proc_acor) 'autocorrelated', if(!err_proc_iid && !err_proc_acor) 'no'), collapse=' and '), ' process error'),
      comment('* ', 'reaeration depends on ',deficit_src),
      
      # process error (always looped)
      if(err_proc_acor) c(
        p(''),
        s('err_proc_acor[1] <- err_proc_acor_inc[1]'),
        p('for(i in 1:(n-2)) {'),
        s('  err_proc_acor[i+1] <- err_proc_acor_phi * err_proc_acor[i] + err_proc_acor_inc[i+1]'),
        p('}')
      ),
      
      # dDO model - applies to any model with deficit_src == 'DO_obs' (includes
      # all process-only models because we've banned 'DO_mod' for them)
      if(deficit_src == 'DO_obs') c(
        p(''),
        comment("dDO model"),
        p('dDO_mod <- '),
        indent(
          if(err_proc_acor) p('err_proc_acor +'),
          p('GPP_daily * coef_GPP +'),
          p('ER_daily * coef_ER +'),
          s('K600_daily * coef_K600_full')
        )
      ),
      
      # DO model - any model that includes observation error
      if(err_obs_iid) c(
        p(''),
        comment("DO model"),
        s('DO_mod[1] <- DO_obs_1'),
        p('for(i in 1:(n-1)) {'),
        indent(
          p('DO_mod[i+1] <- ('),
          p('  DO_mod[i] +'),
          switch(
            deficit_src,
            'DO_obs' = c(
              s('  dDO_mod[i])')),
            'DO_mod' = c(
              p('  GPP_daily * coef_GPP[i] +'),
              p('  ER_daily * coef_ER[i] +'),
              switch(
              ode_method,
              'Euler' = c(
                p('  K600_daily * coef_K600_part[i] * (DO_sat[i] - DO_mod[i])'),
                s(')')),
              'pairmeans' = c(
                p('  K600_daily * coef_K600_part[i] * (DO_sat_pairmean[i] - DO_mod[i]/2)'),
                s(') / (1 + K600_daily * coef_K600_part[i] / 2)'))
              )
            )
          )
        ),
        p('}')
      )
    ),
      
    if(bayes_software == 'stan') c(
      '}',''
    ) else if(bayes_software == 'jags') c(
      ''
    ),
    
    ## model ##
    if(bayes_software == 'stan') c(
      'model {'
    ),
    
    if(err_proc_iid) chunk(
      comment('Independent, identically distributed process error'),
      p('for (i in 1:(n-1)) {'),
      s('  dDO_obs[i] ~ ', f('normal', 'dDO_mod[i]', 'err_proc_iid_sigma')),
      p('}'),
      s('err_proc_iid_sigma ~ ', f('uniform', 'err_proc_iid_sigma_min', 'err_proc_iid_sigma_max'))),
    
    if(err_proc_acor) chunk(
      comment('Autocorrelated process error'),
      p('for(i in 1:(n-1)) {'),
      s('  err_proc_acor_inc[i] ~ ', f('normal', '0', 'err_proc_acor_sigma')),
      p('}'),
      comment('Autocorrelation (phi) & SD (sigma) of the process errors'),
      s('err_proc_acor_phi ~ ', f('uniform', 'err_proc_acor_phi_min', 'err_proc_acor_phi_max')),
      s('err_proc_acor_sigma ~ ', f('uniform', 'err_proc_acor_sigma_min', 'err_proc_acor_sigma_max'))),
    
    if(err_obs_iid) chunk(
      comment('Independent, identically distributed observation error'),
      p('for(i in 1:n) {'),
      s('  DO_obs[i] ~ ', f('normal', 'DO_mod[i]', 'err_obs_iid_sigma')),
      p('}'),
      comment('SD (sigma) of the observation errors'),
      s('err_obs_iid_sigma ~ ', f('uniform', 'err_obs_iid_sigma_min', 'err_obs_iid_sigma_max'))),
    
    indent(
      comment('Daily metabolism values'),
      s('GPP_daily ~ ', f('normal', 'GPP_daily_mu', 'GPP_daily_sigma')),
      s('ER_daily ~ ', f('normal', 'ER_daily_mu', 'ER_daily_sigma')),
      s('K600_daily ~ ', f('normal', 'K600_daily_mu', 'K600_daily_sigma'))),
    
    '}'
  )
  
  writeLines(model_text, con=paste0('inst/models/', model_name), sep="\n")
  
}

#' Generate MCMC code files with all of the desired combinations of features
#' 
#' This function gets run on package build and creates every model within the 
#' set of factorial combinations of arguments to mm_generate_mcmc_file, with the
#' exception of the one pair of incompatible arguments (err_obs_iid=F &&
#' deficit_src='DO_mod')
#' 
#' @keywords internal
mm_generate_mcmc_files <- function() {
  opts <- expand.grid(
    bayes_software=c('jags','stan'), 
    ode_method=c('Euler','pairmeans'),
    deficit_src=c('DO_mod','DO_obs'),
    err_obs_iid=c(TRUE, FALSE),
    err_proc_acor=c(TRUE, FALSE),
    err_proc_iid=c(TRUE, FALSE),
    stringsAsFactors=FALSE)
  attr(opts, 'out.attrs') <- NULL
  
  incompatible <- !opts$err_obs_iid & opts$deficit_src == 'DO_mod'
  opts <- opts[!incompatible, ]
  
  for(i in 1:nrow(opts)) {
    do.call(mm_generate_mcmc_file, opts[i,])
  }
}
mm_generate_mcmc_files()
