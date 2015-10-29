#' Generate the models in inst/models/bayes
#' 
#' @param bayes_software Which software are we generating code for?
#' @param ode_method The method to use in solving the ordinary differential 
#'   equation for DO. Euler: dDOdt from t=1 to t=2 is solely a function of GPP, 
#'   ER, DO, etc. at t=1. pairmeans: dDOdt from t=1 to t=2 is a function of the
#'   mean values of GPP, ER, etc. across t=1 and t=2.
#' @param deficit_src From what DO estimate (observed or modeled) should the DO 
#'   deficit be computed?
#' @keywords internal
mm_generate_mcmc_file <- function(
  bayes_software=c('jags','stan'), 
  ode_method=c('Euler','pairmeans'),
  deficit_src=c('DO_mod','DO_obs'),
  pooling='none') {
  
  # choose/check arguments
  bayes_software <- match.arg(bayes_software)
  ode_method <- match.arg(ode_method)
  deficit_src <- match.arg(deficit_src)
  
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
  
  # name the model
  model_name <- paste0(
    c(none='np', partial='pp')[[pooling]], '_',
    'oi', '_',
    c(Euler='eu', pairmeans='pm')[[ode_method]], '_',
    c(DO_mod='km', DO_obs='ko')[[deficit_src]], '.',
    bayes_software
  )
  
  # define the model
  model_text <- c(
    
    comment(model_name), '',
    
    ## data ##
    if(bayes_software == 'stan') c(
      
      'data {','',
      
      chunk(
        'real GPP_daily_mu;',
        'real GPP_daily_sigma;',
        'real ER_daily_mu;',
        'real ER_daily_sigma;',
        'real K600_daily_mu;',
        'real K600_daily_sigma;'),
      
      chunk(
        'real err_obs_iid_sigma_min;',
        'real err_obs_iid_sigma_max;'),
      
      chunk(
        'int <lower=0> n;',
        'real DO_obs_1;'),
      
      chunk(
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
      'transformed data {',''
    ) else if(bayes_software == 'jags') c(
      'data {',''
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
        if(ode_method == 'pairmeans' && deficit_src == 'DO_mod') 'vector [n-1] DO_sat_pairmean;')
    ),
    
    if(ode_method == 'Euler') { chunk(
      comment('Coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)'),
      p('for(i in 1:(n-1)) {'),
      s('  coef_GPP[i]  <- frac_GPP[i] / depth[i]'),
      s('  coef_ER[i]   <- frac_ER[ i] / depth[i]'),
      switch(
        deficit_src,
        DO_mod = c(
          s('  coef_K600_part[i] <- KO2_conv[i] * frac_D[i]')
        ),
        DO_obs = c(
          p('  coef_K600_full[i] <- KO2_conv[i] * frac_D[i] * '),
          s('    (DO_sat[i] - DO_obs[i])')
        )
      ),
      p('}')
    )} else if(ode_method == 'pairmeans') { chunk(
      comment('Coefficients by pairmeans (e.g., mean(frac_GPP[i:(i+1)]) applies to the DO step from i to i+1)'),
      p('for(i in 1:(n-1)) {'),
      s('  coef_GPP[i]  <- (frac_GPP[i] + frac_GPP[i+1])/2 / ((depth[i] + depth[i+1])/2)'),
      s('  coef_ER[i]   <- (frac_ER[ i] + frac_ER[ i+1])/2 / ((depth[i] + depth[i+1])/2)'),
      switch(
        deficit_src,
        DO_mod = c(
          s('  coef_K600_part[i] <- (KO2_conv[i] + KO2_conv[i+1])/2 * (frac_D[i] + frac_D[i+1])/2'),
          s('  DO_sat_pairmean[i] <- (DO_sat[i] + DO_sat[i+1])/2')
        ),
        DO_obs = c(
          p('  coef_K600_full[i] <- (KO2_conv[i] + KO2_conv[i+1])/2 * (frac_D[i] + frac_D[i+1])/2 *'),
          s('    (DO_sat[i] + DO_sat[i+1] - DO_obs[i] - DO_obs[i+1])/2')
        )
      ),
      p('}')
    )},
    
    # close out transformed data (stan) or data (jags)
    c(
      '}',''
    ),
    
    ## Stan: parameters ##
    if(bayes_software == 'stan') c(
      'parameters {','',
      
      chunk(
        'real GPP_daily;',
        'real ER_daily;',
        'real K600_daily;'), # real DO_mod_1;'
      
      chunk(
        'real <lower=err_obs_iid_sigma_min, upper=err_obs_iid_sigma_max> err_obs_iid_sigma;'),
      
      '}',''
    ),
    
    ## Stan: transformed parameters ## - statements evaluated once per leapfrog step
    ## JAGS: model ##
    if(bayes_software == 'stan') c(
      'transformed parameters {',''
    ) else if(bayes_software == 'jags') c(
      'model {', ''
    ),
    
    if(bayes_software == 'stan') c(
      chunk(
        'vector [n] DO_mod;')
    ),
    
    chunk(
      comment('Model DO time series (',ode_method,' version)'),
      s('DO_mod[1] <- DO_obs_1'), # DO_obs[1] or DO_mod_1;
      p('for(i in 1:(n-1)) {'),
      p('  DO_mod[i+1] <- ('),
      #if(err_proc_acor) p('    err_proc_acor[i] +'),
      p('    DO_mod[i] +'),
      p('    GPP_daily * coef_GPP[i] +'),
      p('    ER_daily * coef_ER[i] +'),
      if(deficit_src == 'DO_obs') { c(
        p('    K600_daily * coef_K600_full[i]'),
        s('  )')
      )} else if(deficit_src == 'DO_mod') {
        if(ode_method == 'Euler') { c(
          p('    K600_daily * coef_K600_part[i] * (DO_sat[i] - DO_mod[i])'),
          s('  )')
        )} else if(ode_method == 'pairmeans') { c(
          p('    K600_daily * coef_K600_part[i] * (DO_sat_pairmean[i] - DO_mod[i]/2)'),
          s('  ) / (1 + K600_daily * coef_K600_part[i] / 2)')
        )}
      },
      p('}')),
    
    if(bayes_software == 'stan') c(
      '}',''
    ),
    
    ## model ##
    if(bayes_software == 'stan') c(
      'model {',''
    ),
    
    chunk(
      comment('Priors on observation error'),
      #s('DO_mod[1] ~ ', f('normal', 'DO_obs_1', 'err_obs_iid_sigma')), # DO_obs[1] or DO_mod_1;
      p('for(i in 1:n) {'),
      s('  DO_obs[i] ~ ', f('normal', 'DO_mod[i]', 'err_obs_iid_sigma')),
      p('}'),
      s('err_obs_iid_sigma ~ ', f('uniform', 'err_obs_iid_sigma_min', 'err_obs_iid_sigma_max'))),
    
    chunk(
      comment('Priors on daily metabolism values'),
      s('GPP_daily ~ ', f('normal', 'GPP_daily_mu', 'GPP_daily_sigma')),
      s('ER_daily ~ ', f('normal', 'ER_daily_mu', 'ER_daily_sigma')),
      s('K600_daily ~ ', f('normal', 'K600_daily_mu', 'K600_daily_sigma'))),
    
    '}'
  )
  
  writeLines(model_text, con=paste0('inst/models/', model_name), sep="\n")
  
}

#' Generate MCMC code files with all of the desired combinations of features
#' 
#' This function gets run on package build
#' 
#' @keywords internal
mm_generate_mcmc_files <- function() {
  opts <- expand.grid(
    bayes_software=c('jags','stan'), 
    ode_method=c('Euler','pairmeans'),
    deficit_src=c('DO_mod','DO_obs'),
    stringsAsFactors=FALSE)
  attr(opts, 'out.attrs') <- NULL
  for(i in 1:length(opts[[1]])) {
    do.call(mm_generate_mcmc_file, opts[i,])
  }
}
mm_generate_mcmc_files()