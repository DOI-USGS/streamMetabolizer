#' Generate the models in inst/models/bayes
#' 
#' @inheritParams mm_name
#' @keywords internal
mm_generate_mcmc_file <- function(
  type='bayes', 
  pool_K600=c('none','normal','linear','binned'),
  err_obs_iid=c(TRUE, FALSE),
  err_proc_acor=c(FALSE, TRUE),
  err_proc_iid=c(FALSE, TRUE),
  ode_method=c('trapezoid','euler'),
  GPP_fun=c('linlight'), #'satlight'
  ER_fun=c('constant'), #'q10temp'
  deficit_src=c('DO_mod','DO_obs'),
  engine='stan') {
  
  # handle Euler and pairmeans as deprecated arguments. mm_name runs a similar check & warning
  if(ode_method %in% c('Euler','pairmeans'))
    warning("for ode_method, 'Euler' and 'pairmeans' are deprecated in favor of 'euler' and 'trapezoid'")
  if(ode_method == 'Euler') ode_method <- 'euler'
  if(ode_method == 'pairmeans') ode_method <- 'trapezoid'
  
  # name the model. much argument checking happens here, even with
  # check_validity=FALSE (which is needed to avoid circularity)
  model_name <- mm_name(
    type='bayes',
    pool_K600=pool_K600,
    err_obs_iid=err_obs_iid, err_proc_acor=err_proc_acor, err_proc_iid=err_proc_iid,
    ode_method=ode_method, GPP_fun=GPP_fun, ER_fun=ER_fun, deficit_src=deficit_src, engine=engine,
    check_validity=FALSE)
  
  # define rules for when to model as process, obs, or both
  dDO_model <- (deficit_src == 'DO_obs')
  DO_model <- (err_obs_iid || deficit_src == 'DO_mod')
  
  #### helper functions ####
  comment <- function(...) { 
    # prefix with the appropriate comment character[s]
    paste0('// ', paste0(...))
  }
  chunk <- function(..., indent=2, newline=TRUE) { 
    # indent a chunk & add a newline
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
  f <- function(distrib, ...) { 
    # check that I've named the arguments correctly
    args <- c(list(...))
    switch(
      distrib,
      normal = { if(!all(names(args) == c('mu','sigma'))) stop("expecting normal(mu,sigma)") }, 
      uniform = { if(!all(names(args) == c('min','max'))) stop("expecting uniform(min,max)") },
      beta = { if(!all(names(args) == c('alpha','beta'))) stop("expecting beta(alpha,beta)") },
      gamma = { if(!all(names(args) == c('shape','rate'))) stop("expecting gamma(shape,rate)")
        # shape = alpha = k = first argument
        # rate = beta = 1/theta = inverse scale = second argument
      },
      lognormal = { if(!all(names(args) == c('location','scale'))) stop("expecting lognormal(location,scale)") 
        # location = mu = first argument
        # scale = sigma = second argument
      }
    )
    # create the function call text
    paste0(distrib, '(', paste0(args, collapse=', '), ')')
  }
  fs <- function(distrib, Y) {
    # Y_scaled is from a standard normal distribution, e.g., Y_scaled ~
    # normal(0, 1). this function returns an equation calculating Y, the
    # rescaled values of Y_scaled that follow the distrib distribution. It
    # assumes there exist parameters with suffixes corresponding to the args
    # required by f()
    paste0(
      Y, ' = ',
      switch(
        distrib,
        normal = stop(),
        uniform = stop(),
        beta = stop(),
        gamma = stop(),
        lognormal = paste0('exp(', paste0(Y, '_location'), ') * pow(exp(', paste0(Y, '_scaled'), '), ', paste0(Y, '_scale'), ')')
      )
    )
  }
  p <- paste0 # partial line: just combine strings into string
  s <- function(...) {
    # line with stop/semicolon: combine strings into string & add semicolon
    p(p(...), ';')
  }
  
  #### <begin model definition> ####
  model_text <- c(
    
    comment(model_name), '',
    
    #### data ####
    c('data {',
      chunk(
        comment('Parameters of priors on metabolism'),
        'real GPP_daily_mu;',
        'real GPP_daily_sigma;',
        'real ER_daily_mu;',
        'real ER_daily_sigma;',
        
        if(pool_K600 %in% c('normal','linear','binned')) c(
          p(''),
          comment('Parameters of hierarchical priors on K600_daily (', pool_K600, ' model)')
        ),
        if(pool_K600 == 'none') c(
          'real K600_daily_mu;',
          'real K600_daily_sigma;'
        ) else c(
          # hierarchical models each do K600_daily_mu/K600_daily_pred
          # differently, but K600_daily_sigma they all do the same
          switch(
            pool_K600,
            normal=c(
              'real K600_daily_mu_mu;',
              'real K600_daily_mu_sigma;'),
            linear=c(
              'vector[2] K600_daily_beta_mu;',
              'vector[2] K600_daily_beta_sigma;'),
            binned=c(
              'int <lower=1> b; # number of K600_daily_betas',
              'vector[b] K600_daily_beta_mu;',
              'vector[b] K600_daily_beta_sigma;')
          ),
          'real K600_daily_sigma_location;',
          'real K600_daily_sigma_scale;'
        )),
      
      chunk(
        comment('Error distributions'),
        if(err_obs_iid) c(
          'real err_obs_iid_sigma_location;',
          'real err_obs_iid_sigma_scale;'),
        if(err_proc_acor) c(
          'real err_proc_acor_phi_alpha;',
          'real err_proc_acor_phi_beta;',
          'real err_proc_acor_sigma_location;',
          'real err_proc_acor_sigma_scale;'),
        if(err_proc_iid) c(
          'real err_proc_iid_sigma_location;',
          'real err_proc_iid_sigma_scale;')),
      
      chunk(
        comment('Data dimensions'),
        'int<lower=1> d; # number of dates',
        'int<lower=1> n; # number of observations per date'),
      
      chunk(
        comment('Daily data'),
        'vector[d] DO_obs_1;',
        switch(
          pool_K600,
          linear='vector[d] ln_discharge_daily;',
          binned='int<lower=1,upper=b> discharge_bin_daily[d];')),
      
      # prepare to iterate over n obs for all d at a time:
      # https://groups.google.com/forum/#!topic/stan-users/ZHeFFV4q_gk
      indent(
        comment('Data'),
        'vector[d] DO_obs[n];',
        'vector[d] DO_sat[n];',
        'vector[d] frac_GPP[n];',
        'vector[d] frac_ER[n];',
        'vector[d] frac_D[n];',
        'vector[d] depth[n];',
        'vector[d] KO2_conv[n];'),
      
      '}',''
    ),
    
    #### transformed data ####
    c('transformed data {', # transformed data = statements evaluated exactly once
      chunk(
        'vector[d] coef_GPP[n-1];',
        'vector[d] coef_ER[n-1];',
        switch(
          deficit_src,
          DO_mod = 'vector[d] coef_K600_part[n-1];',
          DO_obs = 'vector[d] coef_K600_full[n-1];'
        ),
        if(ode_method == 'trapezoid' && deficit_src == 'DO_mod') 'vector[d] DO_sat_pairmean[n-1];',
        if(dDO_model) 'vector[d] dDO_obs[n-1];'
      ),
      
      indent(
        p('for(i in 1:(n-1)) {'),
        indent(
          # Coefficient pre-calculations
          if(ode_method == 'euler') c(
            comment('Coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)'),
            s('coef_GPP[i]  = frac_GPP[i] ./ depth[i]'),
            s('coef_ER[i]   = frac_ER[i] ./ depth[i]'),
            switch(
              deficit_src,
              DO_mod = c(
                s('coef_K600_part[i] = KO2_conv[i] .* frac_D[i]')
              ),
              DO_obs = c(
                p('coef_K600_full[i] = KO2_conv[i] .* frac_D[i] .*'),
                s('  (DO_sat[i] - DO_obs[i])')
              )
            )
          ) else if(ode_method == 'trapezoid') c(
            comment('Coefficients for trapezoid rule (e.g., mean(frac_GPP[i:(i+1)]) applies to the DO step from i to i+1)'),
            s('coef_GPP[i] = (frac_GPP[i] + frac_GPP[i+1])/2.0 ./ ((depth[i] + depth[i+1])/2.0)'),
            s('coef_ER[i] = (frac_ER[i] + frac_ER[i+1])/2.0 ./ ((depth[i] + depth[i+1])/2.0)'),
            switch(
              deficit_src,
              DO_mod = c(
                s('coef_K600_part[i] = (KO2_conv[i] + KO2_conv[i+1])/2.0 .* (frac_D[i] + frac_D[i+1])/2.0'),
                s('DO_sat_pairmean[i] = (DO_sat[i] + DO_sat[i+1])/2.0')
              ),
              DO_obs = c(
                p('coef_K600_full[i] = (KO2_conv[i] + KO2_conv[i+1])/2.0 .* (frac_D[i] + frac_D[i+1])/2.0 .*'),
                s('  (DO_sat[i] + DO_sat[i+1] - DO_obs[i] - DO_obs[i+1])/2.0')
              )
            )
          ),
          
          # dDO pre-calculations
          if(dDO_model) c(
            comment('dDO observations'),
            s('dDO_obs[i] = DO_obs[i+1] - DO_obs[i]')
          )
        ),
        p('}')
      ),
      '}',''
    ),
    
    #### parameters ####
    c('parameters {',
      indent(
        # daily metabolism rate parameters
        c('vector[d] GPP_daily;',
          'vector[d] ER_daily;',
          'vector<lower=0>[d] K600_daily;'),
        
        # K600 pooling parameters
        if(pool_K600 != 'none') c(
          '',
          switch(
            pool_K600,
            normal='real K600_daily_mu;',
            linear='vector[2] K600_daily_beta;',
            binned='vector[b] K600_daily_beta;'),
          'real<lower=0> K600_daily_sigma_scaled;'),
        
        # error distributions
        '',
        if(err_obs_iid) c(
          'real<lower=0> err_obs_iid_sigma_scaled;'),
        if(err_proc_acor) c(
          'real<lower=0, upper=1> err_proc_acor_phi;', # need to figure out how to scale phi (which might be 0-1 or very close to 0)
          'real<lower=0> err_proc_acor_sigma_scaled;'),
        if(err_proc_iid) c(
          'real<lower=0> err_proc_iid_sigma_scaled;'),
        
        # instantaneous process error values
        if((err_proc_iid && !dDO_model) || err_proc_acor) c(
          '',
          if(err_proc_iid && !dDO_model) c(
            'vector[d] err_proc_iid[n-1];'),
          if(err_proc_acor) c(
            'vector[d] err_proc_acor_inc[n-1];'))
      ),
      '}',''
    ),
    
    #### transformed parameters ####
    'transformed parameters {', # transformed parameters = statements evaluated once per leapfrog step
    
    # transformed parameter declarations
    chunk(
      # rescaled K600 pooling parameters
      if(pool_K600 != 'none') c(
        'real K600_daily_sigma;',
        if(pool_K600 %in% c('linear','binned'))
          'vector[d] K600_daily_pred;'
      ),
      
      # rescaled error distribution parameters
      if(err_obs_iid) c(
        'real<lower=0> err_obs_iid_sigma;'),
      if(err_proc_acor) c(
        # 'real<lower=0, upper=1> err_proc_acor_phi;', # need to figure out how to scale phi (which might be 0-1 or very close to 0)
        'real<lower=0> err_proc_acor_sigma;'),
      if(err_proc_iid) c(
        'real<lower=0> err_proc_iid_sigma;'),

      # instantaneous DO, dDO, and/or process error values
      if(DO_model)
        'vector[d] DO_mod[n];',
      if(dDO_model)
        'vector[d] dDO_mod[n-1];',
      if(err_proc_acor)
        'vector[d] err_proc_acor[n-1];'
    ),
    
    chunk(
      comment('Rescale pooling & error distribution parameters'),
      comment('lnN(location,scale) = exp(location)*(exp(N(0,1))^scale)'),
      
      # rescaled K600 pooling parameters
      if(pool_K600 != 'none') c(
        s('K600_daily_sigma = exp(K600_daily_sigma_location) * pow(exp(K600_daily_sigma_scaled), K600_daily_sigma_scale)')
      ),
      
      # rescaled error distribution parameters
      if(err_obs_iid) c(
        s(fs('lognormal', 'err_obs_iid_sigma'))),
      if(err_proc_acor) c(
        # s(fs('beta', 'err_proc_acor_phi'?)), # need to figure out how to scale phi (which might be 0-1 or very close to 0)
        s(fs('lognormal', 'err_proc_acor_sigma'))),
      if(err_proc_iid) c(
        s(fs('lognormal', 'err_proc_iid_sigma')))
    ),
    
    # K600_daily model
    if(pool_K600 %in% c('linear','binned')) chunk(
      comment('Hierarchical, ', pool_K600, ' model of K600_daily'),
      switch(
        pool_K600,
        linear=s('K600_daily_pred = K600_daily_beta[1] + K600_daily_beta[2] * ln_discharge_daily'),
        binned=s('K600_daily_pred = K600_daily_beta[discharge_bin_daily]')
      )
    ),
    
    # model instantaneous DO, dDO, and/or process error values
    indent(
      comment('Model DO time series'),
      comment('* ', ode_method,' version'),
      comment('* ', if(!err_obs_iid) 'no ', 'observation error'),
      comment('* ', paste0(c(if(err_proc_iid) 'IID', if(err_proc_acor) 'autocorrelated', if(!err_proc_iid && !err_proc_acor) 'no'), collapse=' and '), ' process error'),
      comment('* ', 'reaeration depends on ',deficit_src),
      
      # process error (always looped, vectorized across days)
      if(err_proc_acor) c(
        p(''),
        s('err_proc_acor[1] = err_proc_acor_inc[1]'),
        p('for(i in 1:(n-2)) {'),
        s('  err_proc_acor[i+1] = err_proc_acor_phi * err_proc_acor[i] + err_proc_acor_inc[i+1]'),
        p('}')
      ),
      
      # dDO model - applies to any model with deficit_src == 'DO_obs' (includes 
      # all process-only models because we've banned 'DO_mod' for them). looping
      # over times of day, doing all days at once for each time of day
      if(dDO_model) c(
        p(''),
        comment("dDO model"),
        p('for(i in 1:(n-1)) {'),
        indent(
          p('dDO_mod[i] = '),
          indent(
            if(err_proc_acor) p('err_proc_acor[i] +'),
            p('GPP_daily  .* coef_GPP[i] +'),
            p('ER_daily   .* coef_ER[i] +'),
            s('K600_daily .* coef_K600_full[i]')
          )
        ),
        p('}')
      ),
      
      # DO model - any model that includes observation error or is a function of
      # the previous moment's DO_mod
      if(DO_model) c(
        p(''),
        comment("DO model"),
        s('DO_mod[1] = DO_obs_1'),
        p('for(i in 1:(n-1)) {'),
        indent(
          p('DO_mod[i+1] = ('),
          p('  DO_mod[i] +'),
          if(dDO_model) c(
            s('  dDO_mod[i]', if(err_proc_iid) ' + err_proc_iid[i]', ')')
          ) else c(
            if(err_proc_iid) p('  err_proc_iid[i] +'),
            if(err_proc_acor) p('  err_proc_acor[i] +'),
            p('  GPP_daily .* coef_GPP[i] +'),
            p('  ER_daily .* coef_ER[i] +'),
            switch(
              ode_method,
              'euler' = c(
                s('  K600_daily .* coef_K600_part[i] .* (DO_sat[i] - DO_mod[i]))')),
              'trapezoid' = c(
                p('  K600_daily .* coef_K600_part[i] .* (DO_sat_pairmean[i] - DO_mod[i]/2.0)'),
                s(') ./ (1.0 + K600_daily .* coef_K600_part[i] / 2.0)'))
            )
          )
        ),
        p('}')
      )
    ),
    '}','',
    
    #### model ####
    'model {',
    
    if(err_proc_iid || err_proc_acor) chunk(
      comment('Process error'),
      p('for(i in 1:(n-1)) {'),
      indent(
        if(err_proc_iid) c(
          comment('Independent, identically distributed process error'),
          if(dDO_model) s(
            'dDO_obs[i] ~ ', f('normal', mu=p('dDO_mod[i]'), sigma='err_proc_iid_sigma')
          ) else s(
            'err_proc_iid[i] ~ ', f('normal', mu='0', sigma='err_proc_iid_sigma')
          )
        ),
        if(err_proc_acor) c(
          comment('Autocorrelated process error'),
          s('err_proc_acor_inc[i] ~ ', f('normal', mu='0', sigma='err_proc_acor_sigma'))
        )
      ),
      p('}'),
      if(err_proc_iid) c(
        comment('SD (sigma) of the IID process errors'),
        s('err_proc_iid_sigma_scaled ~ ', f('normal', mu='0', sigma='1'))),
      if(err_proc_acor) c(
        comment('Autocorrelation (phi) & SD (sigma) of the process errors'),
        s('err_proc_acor_phi ~ ', f('beta', alpha='err_proc_acor_phi_alpha', beta='err_proc_acor_phi_beta')),
        s('err_proc_acor_sigma_scaled ~ ', f('normal', mu='0', sigma='1')))
    ),
    
    if(err_obs_iid) chunk(
      comment('Independent, identically distributed observation error'),
      #s('DO_mod[1] ~ normal(DO_obs_1, err_obs_iid_sigma)'),
      p('for(i in 2:n) {'),
      indent(
        s('DO_obs[i] ~ ', f('normal', mu=p('DO_mod[i]'), sigma='err_obs_iid_sigma'))
      ),
      p('}'),
      comment('SD (sigma) of the observation errors'),
      s('err_obs_iid_sigma_scaled ~ ', f('normal', mu='0', sigma='1'))),
    
    indent(
      comment('Daily metabolism priors'),
      s('GPP_daily ~ ', f('normal', mu='GPP_daily_mu', sigma='GPP_daily_sigma')),
      s('ER_daily ~ ', f('normal', mu='ER_daily_mu', sigma='ER_daily_sigma')),
      if(pool_K600 %in% c('none','normal')) s(
        'K600_daily ~ ', f('normal', mu='K600_daily_mu', sigma='K600_daily_sigma')
      ) else if(pool_K600 %in% c('linear','binned')) s(
        'K600_daily ~ ', f('normal', mu='K600_daily_pred', sigma='K600_daily_sigma')
      )
    ),
    
    if(pool_K600 != 'none') c(
      '',
      indent(
        comment('Hierarchical constraints on K600_daily (', pool_K600, ' model)'),
        if(pool_K600 == 'normal') c(
          s('K600_daily_mu ~ ', f('normal', mu='K600_daily_mu_mu', sigma='K600_daily_mu_sigma'))
        ),
        if(pool_K600 %in% c('linear','binned')) c(
          s('K600_daily_beta ~ ', f('normal', mu='K600_daily_beta_mu', sigma='K600_daily_beta_sigma'))
        ),
        s('K600_daily_sigma_scaled ~ ', f('normal', mu='0', sigma='1'))
      )
    ),
    
    '}'
  )
  #### <end model definition> ####
  
  writeLines(model_text, con=paste0('inst/models/', model_name), sep="\n")
  
}

#' Generate MCMC code files with all of the desired combinations of features
#' 
#' This function gets run on package build and creates every model within the 
#' set of factorial combinations of arguments to mm_generate_mcmc_file, with the
#' exception of the one pair of incompatible arguments (err_obs_iid=F &&
#' deficit_src='DO_mod')
#' 
#' @include mm_name.R
#' @include mm_parse_name.R
#' @include mm_valid_names.R
#' @include mm_validate_name.R
#' @keywords internal
mm_generate_mcmc_files <- function() {
  opts <- expand.grid(
    pool_K600=c('none','normal','linear','binned'),
    err_obs_iid=c(TRUE, FALSE),
    err_proc_acor=c(FALSE, TRUE),
    err_proc_iid=c(FALSE, TRUE),
    ode_method=c('trapezoid','euler'),
    GPP_fun='linlight',
    ER_fun='constant',
    deficit_src=c('DO_mod','DO_obs'),
    engine='stan',
    stringsAsFactors=FALSE)
  attr(opts, 'out.attrs') <- NULL
  
  incompatible <- 
    (!opts$err_obs_iid & opts$deficit_src == 'DO_mod') |
    (!opts$err_obs_iid & !opts$err_proc_acor & !opts$err_proc_iid)
  opts <- opts[!incompatible, ]
  
  for(i in 1:nrow(opts)) {
    do.call(mm_generate_mcmc_file, opts[i,])
  }
}
mm_generate_mcmc_files()
