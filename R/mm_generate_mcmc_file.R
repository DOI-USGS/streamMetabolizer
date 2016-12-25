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
      beta = { if(!all(names(args) == c('alpha','beta'))) stop("expecting beta(alpha,beta)") },
      gamma = { if(!all(names(args) == c('shape','rate'))) stop("expecting gamma(shape,rate)")
        # shape = alpha = k = first argument
        # rate = beta = 1/theta = inverse scale = second argument
      },
      halfcauchy = { if(!all(names(args) == c('scale'))) stop("expecting halfcauchy(scale)")
        distrib <- 'cauchy'
        args <- c(list(location=0), args)
      },
      halfnormal = { if(!all(names(args) == c('sigma'))) stop("expecting halfnormal(sigma)")
        distrib <- 'normal'
        args <- c(list(mu=0), args)
      },
      lognormal = { if(!all(names(args) == c('meanlog','sdlog'))) stop("expecting lognormal(meanlog,sdlog)") 
        # meanlog = mu = first argument
        # sdlog = sigma = second argument
      },
      normal = { if(!all(names(args) == c('mu','sigma'))) stop("expecting normal(mu,sigma)") }, 
      uniform = { if(!all(names(args) == c('min','max'))) stop("expecting uniform(min,max)") }
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
        beta = stop(),
        gamma = stop(),
        halfcauchy = sprintf('%s_scale * %s_scaled', Y, Y), # scaled = cauchy(0,1)
        lognormal = sprintf('exp(%s_meanlog + %s_sdlog * %s_scaled)', Y, Y, Y), # scaled = norm(0,1)
        normal = sprintf('%s_sigma * %s_scaled', Y, Y), # scaled = norm(0,1)
        uniform = sprintf('%s_min + (%s_max - %s_min) * %s_scaled', Y, Y) # scaled = unif(0,1)
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
        'real<lower=0> GPP_daily_sigma;',
        'real ER_daily_mu;',
        'real<lower=0> ER_daily_sigma;',
        
        if(pool_K600 %in% c('normal','linear','binned')) c(
          p(''),
          comment('Parameters of hierarchical priors on K600_daily (', pool_K600, ' model)')
        ),
        if(pool_K600 == 'none') c(
          'real K600_daily_meanlog;'
        ) else c(
          # hierarchical models each do K600_daily_meanlog/K600_daily_predlog
          # differently, but K600_daily_sdlog they all do the same
          switch(
            pool_K600,
            normal=c(
              'real K600_daily_meanlog_meanlog;',
              'real<lower=0> K600_daily_meanlog_sdlog;'),
            linear=c(
              'real lnK600_lnQ_intercept_mu;',
              'real<lower=0> lnK600_lnQ_intercept_sigma;',
              'real lnK600_lnQ_slope_mu;',
              'real<lower=0> lnK600_lnQ_slope_sigma;'),
            binned=c(
              'int <lower=1> b; # number of K600_lnQ_nodes',
              'real K600_lnQ_nodediffs_sdlog;',
              'vector[b] K600_lnQ_nodes_meanlog;',
              'vector[b] K600_lnQ_nodes_sdlog;')
          )
        ),
        'real<lower=0> K600_daily_sdlog;'
        ),
      
      chunk(
        comment('Error distributions'),
        if(err_obs_iid) c(
          'real<lower=0> err_obs_iid_sigma_scale;'),
        if(err_proc_acor) c(
          'real err_proc_acor_phi_alpha;',
          'real err_proc_acor_phi_beta;',
          'real<lower=0> err_proc_acor_sigma_scale;'),
        if(err_proc_iid) c(
          'real<lower=0> err_proc_iid_sigma_scale;')),
      
      chunk(
        comment('Data dimensions'),
        'int<lower=1> d; # number of dates',
        'real<lower=0> timestep; # length of each timestep in days',
        'int<lower=1> n24; # number of observations in first 24 hours per date',
        'int<lower=1> n; # number of observations per date'),
      
      chunk(
        comment('Daily data'),
        'vector[d] DO_obs_1;',
        switch(
          pool_K600,
          linear=c(
            'vector[d] lnQ_daily;'),
          binned=c(
            'int<lower=1,upper=b> lnQ_bins[2,d];',
            'vector<lower=0,upper=1>[d] lnQ_bin_weights[2];')
        )),
      
      # prepare to iterate over n obs for all d at a time:
      # https://groups.google.com/forum/#!topic/stan-users/ZHeFFV4q_gk
      indent(
        comment('Data'),
        'vector[d] DO_obs[n];',
        'vector[d] DO_sat[n];',
        # as of 10/13/2016 frac_GPP and frac_ER need to be multipliers rather
        # than fractions (i.e., must yield per-day rather than per-timestep
        # rates)
        'vector[d] frac_GPP[n];', 
        'vector[d] frac_ER[n];',
        'vector[d] frac_D[n];',
        'vector[d] depth[n];',
        'vector[d] KO2_conv[n];'),
      
      '}',''
    ),
    
    #### transformed data ####
    # c('transformed data {', # transformed data = statements evaluated exactly once
    #   indent(
    #     #   chunk(
    #     #     # Coefficient declarations, if any, go here
    #     #   ),
    #     #   
    #     #   indent(
    #     #     p('for(i in 1:n) {'),
    #     #     indent(
    #     #       # Coefficient pre-calculations, if any, go here
    #     #     ),
    #     #     p('}')
    #     #   ),
    #   ),
    #   '}',''
    # ),
    
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
            normal=c(
              'real K600_daily_predlog;'),
            linear=c(
              'real lnK600_lnQ_intercept;',
              'real lnK600_lnQ_slope;'),
            binned=c(
              'vector[b] lnK600_lnQ_nodes;')
          )),
        
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
        if(err_proc_acor) c(
          '',
          sprintf('vector[d] err_proc_acor_inc[%s];', switch(ode_method, euler='n-1', trapezoid='n'))),
        
        # DO_mod if it's a fitted parameter (oipi models)
        if(err_obs_iid && err_proc_iid) c(
          'vector[d] DO_mod[n];')
        
      ),
      '}',''
    ),
    
    #### transformed parameters ####
    'transformed parameters {', # transformed parameters = statements evaluated once per leapfrog step
    
    # transformed parameter declarations
    chunk(
      # rescaled K600 pooling parameters
      if(pool_K600 %in% c('linear','binned')) c(
        'vector[d] K600_daily_predlog;'
      ),
      
      # rescaled error distribution parameters
      if(err_obs_iid) c(
        'real<lower=0> err_obs_iid_sigma;'),
      if(err_proc_acor || err_proc_iid) c(
        'vector[d] DO_mod_partial_sigma[n];'
      ),
      if(err_proc_acor) c(
        # 'real<lower=0, upper=1> err_proc_acor_phi;', # currently opting not to scale phi (which might be 0-1 or very close to 0)
        'real<lower=0> err_proc_acor_sigma;'),
      if(err_proc_iid) c(
        'real<lower=0> err_proc_iid_sigma;'),

      # instantaneous GPP, ER, and KO2. the nth value isn't used to calculate DO
      # when ode_method=euler, but it's always used to calculate GPP and ER
      c('vector[d] GPP_inst[n];',
        'vector[d] ER_inst[n];',
        'vector[d] KO2_inst[n];'),
      
      # instantaneous DO and possibly process error values
      if(err_proc_iid)
        'vector[d] DO_mod_partial[n];'
      else # err_obs_iid and/or err_proc_acor without err_proc_iid
        'vector[d] DO_mod[n];',
      if(err_proc_acor)
        sprintf('vector[d] err_proc_acor[%s];', switch(ode_method, euler='n-1', trapezoid='n'))
    ),
    
    # error distribution parameters
    chunk(
      comment('Rescale error distribution parameters'),
      if(err_obs_iid) c(
        s(fs('halfcauchy', 'err_obs_iid_sigma'))),
      if(err_proc_acor) c(
        # s(fs('beta', 'err_proc_acor_phi'?)), # currently opting not to scale phi (which might be 0-1 or very close to 0)
        s(fs('halfcauchy', 'err_proc_acor_sigma'))),
      if(err_proc_iid) c(
        s(fs('halfcauchy', 'err_proc_iid_sigma')))
    ),
    
    # K600_daily model
    if(pool_K600 %in% c('linear','binned')) chunk(
      comment('Hierarchical, ', pool_K600, ' model of K600_daily'),
      switch(
        pool_K600,
        linear=s('K600_daily_predlog = lnK600_lnQ_intercept + lnK600_lnQ_slope * lnQ_daily'),
        binned=s('K600_daily_predlog = lnK600_lnQ_nodes[lnQ_bins[1]] .* lnQ_bin_weights[1] + \n  ',
                 '                     lnK600_lnQ_nodes[lnQ_bins[2]] .* lnQ_bin_weights[2]')
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
        comment("Calculate autocorrelated process error rates"),
        s('err_proc_acor[1] = err_proc_acor_inc[1]'),
        p(sprintf('for(i in 2:%s) {', switch(ode_method, euler='(n-1)', trapezoid='n'))),
        s('  err_proc_acor[i] = err_proc_acor_phi * err_proc_acor[i-1] + err_proc_acor_inc[i]'),
        p('}')
      ),

      # individual processes
      c(
        p(''),
        comment("Calculate individual process rates"),
        p('for(i in 1:n) {'),
        indent(
          s('GPP_inst[i] = GPP_daily .* frac_GPP[i]'),
          s('ER_inst[i] = ER_daily .* frac_ER[i]'),
          s('KO2_inst[i] = K600_daily .* KO2_conv[i]')
        ),
        p('}')
      ),
      
      # DO model - any model that includes observation error or is a function of
      # the previous moment's DO_mod
      c(
        p(''),
        comment("DO model"),
        if(!err_proc_iid)
          s('DO_mod[1] = DO_obs_1')
        else if(err_proc_iid && !err_obs_iid && deficit_src=='DO_mod') c(
          # pi_km and pcpi_km models are the only ones where we need
          # DO_mod_partial[1]; others start with [2] or don't use DO_mod_partial
          s('DO_mod_partial[1] = DO_obs_1'),
          s('DO_mod_partial_sigma[1] = err_proc_iid_sigma * timestep ./ depth[1]')
        ),
        p('for(i in 1:(n-1)) {'),
        indent(
          if(err_proc_iid) p(
            'DO_mod_partial[i+1] ='
          ) else p( # err_obs_iid and/or err_proc_acor without err_proc_iid
            'DO_mod[i+1] ='
          ),
          indent(
            switch(
              ode_method,
              'euler' = c(
                p(if(err_obs_iid) 'DO_mod' else 'DO_obs', '[i] + (')
              ),
              'trapezoid' = c(
                switch(
                  deficit_src,
                  'DO_obs' = c(
                    p(if(err_obs_iid) 'DO_mod' else 'DO_obs', '[i] + ('),
                    p('  - KO2_inst[i] .* DO_obs[i] - KO2_inst[i+1] .* DO_obs[i+1] +')),
                  'DO_mod' = c(
                    p(if(err_obs_iid) 'DO_mod' else 'DO_mod_partial', '[i] .*'),
                    p('  (2.0 - KO2_inst[i] * timestep) ./ (2.0 + KO2_inst[i+1] * timestep) + ('))
                )
              )
            ),
            p('  (GPP_inst[i] + ER_inst[i]', if(err_proc_acor) ' + err_proc_acor[i]', ') ./ depth[i] +'),
            switch(
              ode_method,
              'euler' = c(
                p('  KO2_inst[i] .* (DO_sat[i] - ',
                  if(deficit_src=='DO_mod' && !err_obs_iid) 'DO_mod_partial' else deficit_src,
                  '[i])'),
                s(') * timestep')
              ),
              'trapezoid' = c(
                p('  (GPP_inst[i+1] + ER_inst[i+1]', if(err_proc_acor) ' + err_proc_acor[i+1]', ') ./ depth[i+1] +'),
                p('  KO2_inst[i] .* DO_sat[i] + KO2_inst[i+1] .* DO_sat[i+1]'),
                switch(
                  deficit_src,
                  'DO_obs' = c(
                    s(') * (timestep / 2.0)')),
                  'DO_mod' = c(
                    s(') .* (timestep ./ (2.0 + KO2_inst[i+1] * timestep))'))
                )
              )
            )
          ),
          if(err_proc_iid) c(
            p('for(j in 1:d) {'),
            indent(
              'DO_mod_partial_sigma[i+1,j] = err_proc_iid_sigma * ',
              switch(
                ode_method,
                'euler' = indent(
                  s('timestep ./ depth[i,j]')
                ),
                'trapezoid' = indent(
                  'sqrt(pow(depth[i,j], -2) + pow(depth[i+1,j], -2)) .*',
                  switch(
                    deficit_src,
                    'DO_obs' = s('(timestep / 2.0)'),
                    'DO_mod' = s('(timestep / (2.0 + KO2_inst[i+1,j] * timestep))')
                  )
                )
              )
            ),
            p('}')
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
      p('for(i in 2:n) {'),
      indent(
        if(err_proc_iid) c(
          comment('Independent, identically distributed process error'),
          s(if(!err_obs_iid) 'DO_obs[i]' else 'DO_mod[i]', ' ~ ', 
            f('normal', mu='DO_mod_partial[i]', sigma='DO_mod_partial_sigma[i]')
          )
        ),
        if(err_proc_acor) c(
          comment('Autocorrelated process error'),
          s('err_proc_acor_inc[i-1] ~ ', f('normal', mu='0', sigma='err_proc_acor_sigma'))
        )
      ),
      p('}'),
      if(err_proc_iid) c(
        comment('SD (sigma) of the IID process errors'),
        s('err_proc_iid_sigma_scaled ~ ', f('halfcauchy', scale='1'))),
      if(err_proc_acor) c(
        comment('Autocorrelation (phi) & SD (sigma) of the process errors'),
        s('err_proc_acor_phi ~ ', f('beta', alpha='err_proc_acor_phi_alpha', beta='err_proc_acor_phi_beta')), # currently opting not to scale phi (which might be 0-1 or very close to 0)
        s('err_proc_acor_sigma_scaled ~ ', f('halfcauchy', scale='1')))
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
      s('err_obs_iid_sigma_scaled ~ ', f('halfcauchy', scale='1'))),
    
    indent(
      comment('Daily metabolism priors'),
      s('GPP_daily ~ ', f('normal', mu='GPP_daily_mu', sigma='GPP_daily_sigma')),
      s('ER_daily ~ ', f('normal', mu='ER_daily_mu', sigma='ER_daily_sigma')),
      if(pool_K600 %in% c('none')) s(
        'K600_daily ~ ', f('lognormal', meanlog='K600_daily_meanlog', sdlog='K600_daily_sdlog')
      ) else if(pool_K600 %in% c('normal','linear','binned')) s(
        'K600_daily ~ ', f('lognormal', meanlog='K600_daily_predlog', sdlog='K600_daily_sdlog')
      )
    ),
    
    if(pool_K600 != 'none') chunk(
      comment('Hierarchical constraints on K600_daily (', pool_K600, ' model)'),
      switch(
        pool_K600,
        'normal' = c(
          s('K600_daily_predlog ~ ', f('normal', mu='K600_daily_meanlog_meanlog', sigma='K600_daily_meanlog_sdlog '))
        ),
        'linear' = c(
          s('lnK600_lnQ_intercept ~ ', f('normal', mu='lnK600_lnQ_intercept_mu', sigma='lnK600_lnQ_intercept_sigma')),
          s('lnK600_lnQ_slope ~ ', f('normal', mu='lnK600_lnQ_slope_mu', sigma='lnK600_lnQ_slope_sigma'))
        ),
        'binned' = c(
          s('lnK600_lnQ_nodes ~ ', f('normal', mu='K600_lnQ_nodes_meanlog', sigma='K600_lnQ_nodes_sdlog')),
          p('for(k in 2:b) {'),
          s('  lnK600_lnQ_nodes[k] ~ ', f('normal', mu='lnK600_lnQ_nodes[k-1]', sigma='K600_lnQ_nodediffs_sdlog')),
          p('}')
        )
      )
    ),
    
    '}',
    
    #### generated quantities ####
    'generated quantities {',
    
    chunk(
      if(err_obs_iid) 'vector[d] err_obs_iid[n-1];',
      if(err_proc_iid) 'vector[d] err_proc_iid[n-1];',
      'vector[d] GPP;',
      'vector[d] ER;',
      '',
      if(err_obs_iid || err_proc_iid) c(
        'for(i in 1:(n-1)) {',
        indent(
          if(err_obs_iid) 
            s('err_obs_iid[i] = DO_mod[i+1] - DO_obs[i+1]'),
          if(err_proc_iid) 
            s('err_proc_iid[i] = (DO_mod_partial[i+1] - ', if(!err_obs_iid) 'DO_obs[i+1]' else 'DO_mod[i+1]', 
              ') .* (err_proc_iid_sigma ./ DO_mod_partial_sigma[i+1])')
        ),
        p('}')
      ),
      'for(j in 1:d) {',
      indent(
        s('GPP[j] = sum(GPP_inst[1:n24,j]) / n24'),
        s('ER[j] = sum(ER_inst[1:n24,j]) / n24')
      ),
      p('}')
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
    err_proc_acor=FALSE,
    err_proc_iid=c(FALSE, TRUE),
    ode_method=c('trapezoid','euler'),
    GPP_fun='linlight',
    ER_fun='constant',
    deficit_src=c('DO_mod','DO_obs'),
    engine='stan',
    stringsAsFactors=FALSE)
  attr(opts, 'out.attrs') <- NULL
  
  incompatible <- 
    (!opts$err_obs_iid & !opts$err_proc_acor & !opts$err_proc_iid) # need at least 1 kind of error
  opts <- opts[!incompatible, ]
  
  for(i in 1:nrow(opts)) {
    do.call(mm_generate_mcmc_file, opts[i,])
  }
}
mm_generate_mcmc_files()
