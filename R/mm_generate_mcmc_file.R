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
  ode_method=c('pairmeans','Euler'),
  deficit_src=c('DO_mod','DO_obs'),
  engine=c('stan','jags')) {
  
  # name the model. much argument checking happens here, even with
  # check_validity=FALSE (which is needed to avoid circularity)
  model_name <- mm_name(
    type='bayes',
    pool_K600=pool_K600,
    err_obs_iid=err_obs_iid, err_proc_acor=err_proc_acor, err_proc_iid=err_proc_iid,
    ode_method=ode_method, deficit_src=deficit_src, engine=engine,
    check_validity=FALSE)
  
  # define rules for when to model as process, obs, or both
  dDO_model <- (deficit_src == 'DO_obs')
  DO_model <- (err_obs_iid || deficit_src == 'DO_mod')
  
  # helper functions
  comment <- function(...) { 
    # prefix with the appropriate comment character[s]
    chr <- switch(engine, jags='#', stan='//')
    paste0(chr, ' ', paste0(...))
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
  condindent <- function(..., engines=c('jags'), indent=2) {
    # conditional indent; only indent if engine is in engines
    chunk(..., indent=if(engine %in% engines) indent else 0, newline=FALSE)
  }
  condloop <- function(..., iter='j', over='1:d', engines=c('jags'), indent=2) {
    # conditional loop
    do_indent <- engine %in% engines
    c(
      if(do_indent) p('for(',iter,' in ',over,') {'),
      chunk(..., indent=if(do_indent) indent else 0, newline=FALSE),
      if(do_indent) p('}')
    )
  }
  f <- function(distrib, ...) { 
    # switch things to JAGS format if needed
    args <- c(list(...))
    if(engine=='jags') {
      switch(
        distrib,
        normal = {
          if(!all(names(args) == c('mu','sigma'))) stop("expecting normal(mu,sigma)")
          distrib <- 'dnorm'
          args['sigma'] <- sprintf('pow(%s, -2)', args['sigma'])
        }, 
        uniform = {
          if(!all(names(args) == c('min','max'))) stop("expecting uniform(min,max)")
          distrib <- 'dunif'
        },
        gamma = {
          # shape = alpha = k = first argument to both stan and jags
          # rate = beta = 1/theta = inverse scale = second argument to both stan and jags
          if(!all(names(args) == c('shape','rate'))) stop("expecting gamma(shape,rate)")
          distrib <- 'dgamma'
        }
      )
    }
    paste0(
      distrib, '(',
      paste0(args, collapse=', '),
      ')')
  }
  e <- function(operator) {
    # element-wise multiplication or division is .* or ./ in Stan but just * or
    # . in JAGS
    paste0(
      switch(engine, jags='', stan='.'),
      operator)
  }
  p <- paste0 # partial line: just combine strings into string
  s <- function(...) {
    # line with stop/semicolon: combine strings into string & add semicolon if
    # needed
    p(p(...), switch(engine, jags='', stan=';'))
  }
  N <- function(n) {
    # index into a matrix according to the index of the time of day (the col for
    # Stan or the row for JAGS)
    switch(
      engine,
      jags=paste0('[1:d,',n,']'),
      stan=paste0('[',n,']'))
  }
  MN <- function(m, n) {
    # matrix indexing for distributions, which are vectorizable in Stan (go by
    # col) but not in JAGS (go by cell)
    switch(
      engine,
      jags=paste0('[',m,',',n,']'),
      stan=paste0('[',n,']'))
  }
  M <- function(m) {
    # day-by-day or beta-by-beta indexing for distributions of daily values, 
    # which are vectorizable in Stan but not JAGS
    switch(
      engine,
      jags=paste0('[',m,']'),
      stan='')
  }
  
  # define the model
  model_text <- c(
    
    comment(model_name), '',
    
    ## data ##
    if(engine == 'stan') c(
      
      'data {',
      
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
          'real K600_daily_sigma_shape;',
          'real K600_daily_sigma_rate;'
        )),
      
      chunk(
        comment('Error distributions'),
        if(err_obs_iid) c(
          'real err_obs_iid_sigma_shape;',
          'real err_obs_iid_sigma_rate;'),
        if(err_proc_acor) c(
          'real err_proc_acor_phi_shape;',
          'real err_proc_acor_phi_rate;',
          'real err_proc_acor_sigma_shape;',
          'real err_proc_acor_sigma_rate;'),
        if(err_proc_iid) c(
          'real err_proc_iid_sigma_shape;',
          'real err_proc_iid_sigma_rate;')),
      
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
    
    ## Stan: transformed data ## - statements evaluated exactly once
    ## JAGS: data ##
    if(engine == 'stan') c(
      'transformed data {'
    ) else if(engine == 'jags') c(
      'data {'
    ),
    
    if(engine == 'stan') c(
      chunk(
        'vector[d] coef_GPP[n-1];',
        'vector[d] coef_ER[n-1];',
        switch(
          deficit_src,
          DO_mod = 'vector[d] coef_K600_part[n-1];',
          DO_obs = 'vector[d] coef_K600_full[n-1];'
        ),
        if(ode_method == 'pairmeans' && deficit_src == 'DO_mod') 'vector[d] DO_sat_pairmean[n-1];',
        if(dDO_model) 'vector[d] dDO_obs[n-1];')
    ),
    
    indent(
      p('for(i in 1:(n-1)) {'),
      indent(
        # Coefficient pre-calculations
        if(ode_method == 'Euler') c(
          comment('Coefficients by lag (e.g., frac_GPP[i] applies to the DO step from i to i+1)'),
          s('coef_GPP', N('i'), '  <- frac_GPP', N('i'), ' ', e('/'), ' depth', N('i'), ''),
          s('coef_ER', N('i'), '   <- frac_ER',  N('i'), ' ', e('/'), ' depth', N('i'), ''),
          switch(
            deficit_src,
            DO_mod = c(
              s('coef_K600_part', N('i'), ' <- KO2_conv', N('i'), ' ', e('*'), ' frac_D', N('i'), '')
            ),
            DO_obs = c(
              p('coef_K600_full', N('i'), ' <- KO2_conv', N('i'), ' ', e('*'), ' frac_D', N('i'), ' ', e('*')),
              s('  (DO_sat', N('i'), ' - DO_obs', N('i'), ')')
            )
          )
        ) else if(ode_method == 'pairmeans') c(
          comment('Coefficients by pairmeans (e.g., mean(frac_GPP[i:(i+1)]) applies to the DO step from i to i+1)'),
          s('coef_GPP', N('i'), '  <- (frac_GPP', N('i'), ' + frac_GPP', N('i+1'), ')/2.0 ', e('/'), ' ((depth', N('i'), ' + depth', N('i+1'), ')/2.0)'),
          s('coef_ER', N('i'), '   <- (frac_ER',  N('i'), ' + frac_ER',  N('i+1'), ')/2.0 ', e('/'), ' ((depth', N('i'), ' + depth', N('i+1'), ')/2.0)'),
          switch(
            deficit_src,
            DO_mod = c(
              s('coef_K600_part', N('i'), ' <- (KO2_conv', N('i'), ' + KO2_conv', N('i+1'), ')/2.0 ', e('*'), ' (frac_D', N('i'), ' + frac_D', N('i+1'), ')/2.0'),
              s('DO_sat_pairmean', N('i'), ' <- (DO_sat', N('i'), ' + DO_sat', N('i+1'), ')/2.0')
            ),
            DO_obs = c(
              p('coef_K600_full', N('i'), ' <- (KO2_conv', N('i'), ' + KO2_conv', N('i+1'), ')/2.0 ', e('*'), ' (frac_D', N('i'), ' + frac_D', N('i+1'), ')/2.0 ', e('*')),
              s('  (DO_sat', N('i'), ' + DO_sat', N('i+1'), ' - DO_obs', N('i'), ' - DO_obs', N('i+1'), ')/2.0')
            )
          )
        ),
        
        # dDO pre-calculations
        if(dDO_model) c(
          comment('dDO observations'),
          s('dDO_obs', N('i'), ' <- DO_obs', N('i+1'), ' - DO_obs', N('i'), '')
        )
        
        # # vector of ones pre-calculations
        # if(dDO_model && engine == 'jags') c(
        #   comment('ones vector for expanding daily vectors into matrices'),
        #   p('for (i in 1:(n-1)) {'),
        #   indent(s('ones[i] <- 1')),
        #   p('}')
        # )
      ),
      p('}')
    ),
    
    # close out transformed data (stan) or data (jags)
    c(
      '}',''
    ),
    
    ## Stan: parameters ##
    if(engine == 'stan') c(
      'parameters {',
      indent(
        c('vector[d] GPP_daily;',
          'vector[d] ER_daily;',
          'vector[d] K600_daily;'),
        
        if(pool_K600 != 'none') c(
          '',
          switch(
            pool_K600,
            normal='real K600_daily_mu;',
            linear='vector[2] K600_daily_beta;',
            binned='vector[b] K600_daily_beta;'),
          'real K600_daily_sigma;'),
        
        if((err_proc_iid && !dDO_model) || err_proc_acor) c(
          ''),
        if(err_proc_iid && !dDO_model) c(
          'vector[d] err_proc_iid[n-1];'),
        if(err_proc_acor) c(
          'vector[d] err_proc_acor_inc[n-1];'),
        
        if(err_obs_iid || err_proc_acor || err_proc_iid) c(
          '',
          if(err_obs_iid) c(
            'real err_obs_iid_sigma;'),
          if(err_proc_acor) c(
            'real err_proc_acor_phi;',
            'real err_proc_acor_sigma;'),
          if(err_proc_iid) c(
            'real err_proc_iid_sigma;'))
      ),
      '}',''
    ),
    
    ## Stan: transformed parameters ## - statements evaluated once per leapfrog step
    ## JAGS: model ##
    if(engine == 'stan') c(
      'transformed parameters {'
    ) else if(engine == 'jags') c(
      'model {'
    ),
    
    if(engine == 'stan') chunk(
      if(DO_model)
        'vector[d] DO_mod[n];',
      if(dDO_model)
        'vector[d] dDO_mod[n-1];',
      if(err_proc_acor) 
        'vector[d] err_proc_acor[n-1];',
      if(pool_K600 %in% c('linear','binned'))
        'vector[d] K600_daily_pred;'),
    
    indent(
      comment('Model DO time series'),
      comment('* ', ode_method,' version'),
      comment('* ', if(!err_obs_iid) 'no ', 'observation error'),
      comment('* ', paste0(c(if(err_proc_iid) 'IID', if(err_proc_acor) 'autocorrelated', if(!err_proc_iid && !err_proc_acor) 'no'), collapse=' and '), ' process error'),
      comment('* ', 'reaeration depends on ',deficit_src),
      
      # process error (always looped, vectorized across days)
      if(err_proc_acor) c(
        p(''),
        s('err_proc_acor', N('1'), ' <- err_proc_acor_inc', N('1'), ''),
        p('for(i in 1:(n-2)) {'),
        s('  err_proc_acor', N('i+1'), ' <- err_proc_acor_phi * err_proc_acor', N('i'), ' + err_proc_acor_inc', N('i+1'), ''),
        p('}')
      ),
      
      # dDO model - applies to any model with deficit_src == 'DO_obs' (includes 
      # all process-only models because we've banned 'DO_mod' for them). this 
      # should be doable with a single set of elementwise matrix multiplications
      # rather than a loop, but there's no vector[] .* vector[] and no 
      # to_matrix(vector[]), which is slowing me down in Stan. haven't even
      # tried in JAGS.
      if(dDO_model) c(
        p(''),
        comment("dDO model"),
        p('for(i in 1:(n-1)) {'),
        indent(
          p('dDO_mod', N('i'), ' <- '),
          indent(
            if(err_proc_acor) p('err_proc_acor', N('i'), ' +'),
            p('GPP_daily  ', e('*'), ' coef_GPP', N('i'), ' +'),
            p('ER_daily   ', e('*'), ' coef_ER', N('i'), ' +'),
            s('K600_daily ', e('*'), ' coef_K600_full', N('i'), '')
          )
        ),
        p('}')
      ),
      
      # DO model - any model that includes observation error or is a function of
      # the previous moment's DO_mod
      if(DO_model) c(
        p(''),
        comment("DO model"),
        s('DO_mod', N('1'), ' <- DO_obs_1'),
        p('for(i in 1:(n-1)) {'),
        indent(
          p('DO_mod', N('i+1'), ' <- ('),
          p('  DO_mod', N('i'), ' +'),
          if(dDO_model) c(
            s('  dDO_mod', N('i'), ')')
          ) else c(
            if(err_proc_iid) p('  err_proc_iid', N('i'), ' +'),
            if(err_proc_acor) p('  err_proc_acor', N('i'), ' +'),
            p('  GPP_daily ', e('*'), ' coef_GPP', N('i'), ' +'),
            p('  ER_daily ', e('*'), ' coef_ER', N('i'), ' +'),
            switch(
              ode_method,
              'Euler' = c(
                p('  K600_daily ', e('*'), ' coef_K600_part', N('i'), ' ', e('*'), ' (DO_sat', N('i'), ' - DO_mod', N('i'), ')'),
                s(')')),
              'pairmeans' = c(
                p('  K600_daily ', e('*'), ' coef_K600_part', N('i'), ' ', e('*'), ' (DO_sat_pairmean', N('i'), ' - DO_mod', N('i'), '/2.0)'),
                s(') ', e('/'), ' (1.0 + K600_daily ', e('*'), ' coef_K600_part', N('i'), ' / 2.0)'))
            )
          )
        ),
        p('}')
      ),
      
      # K600_daily model
      if(pool_K600 %in% c('linear','binned')) c(
        p(''),
        comment('Hierarchical, ', pool_K600, ' model of K600_daily'),
        switch(
          pool_K600,
          linear=s('K600_daily_pred <- K600_daily_beta[1] + K600_daily_beta[2] * ln_discharge_daily'),
          binned=c(
            # JAGS might not require a loop here
            condloop(
              s('K600_daily_pred', M('j'), ' <- K600_daily_beta[discharge_bin_daily', M('j'), ']')
            )
          )
        )
      )
      
    ),
    
    if(engine == 'stan') c(
      '}',''
    ) else if(engine == 'jags') c(
      ''
    ),
    
    ## Stan: model ##
    if(engine == 'stan') c(
      'model {'
    ),
    
    if(err_proc_iid) chunk(
      comment('Independent, identically distributed process error'),
      p('for(i in 1:(n-1)) {'),
      indent(
        condloop(
          if(dDO_model) s(
            'dDO_obs', MN('j','i'), ' ~ ', f('normal', mu=p('dDO_mod', MN('j','i')), sigma='err_proc_iid_sigma')
          ) else s(
            'err_proc_iid', MN('j','i'), ' ~ ', f('normal', mu='0', sigma='err_proc_iid_sigma')
          )
        )
      ),
      p('}'),
      comment('SD (sigma) of the IID process errors'),
      s('err_proc_iid_sigma ~ ', f('gamma', shape='err_proc_iid_sigma_shape', rate='err_proc_iid_sigma_rate'))),
    
    if(err_proc_acor) chunk(
      comment('Autocorrelated process error'),
      p('for(i in 1:(n-1)) {'),
      indent(
        condloop(
          s('err_proc_acor_inc', MN('j','i'), ' ~ ', f('normal', mu='0', sigma='err_proc_acor_sigma'))
        )
      ),
      p('}'),
      comment('Autocorrelation (phi) & SD (sigma) of the process errors'),
      s('err_proc_acor_phi ~ ', f('gamma', shape='err_proc_acor_phi_shape', rate='err_proc_acor_phi_rate')),
      s('err_proc_acor_sigma ~ ', f('gamma', shape='err_proc_acor_sigma_shape', rate='err_proc_acor_sigma_rate'))),
    
    if(err_obs_iid) chunk(
      comment('Independent, identically distributed observation error'),
      p('for(i in 1:n) {'),
      indent(
        condloop(
          s('DO_obs', MN('j','i'), ' ~ ', f('normal', mu=p('DO_mod', MN('j','i')), sigma='err_obs_iid_sigma'))
        )
      ),
      p('}'),
      comment('SD (sigma) of the observation errors'),
      s('err_obs_iid_sigma ~ ', f('gamma', shape='err_obs_iid_sigma_shape', rate='err_obs_iid_sigma_rate'))),
    
    indent(
      comment('Daily metabolism priors'),
      condloop(
        s('GPP_daily', M('j'), ' ~ ', f('normal', mu='GPP_daily_mu', sigma='GPP_daily_sigma')),
        s('ER_daily', M('j'), ' ~ ', f('normal', mu='ER_daily_mu', sigma='ER_daily_sigma')),
        if(pool_K600 %in% c('none','normal')) s(
          'K600_daily', M('j'), ' ~ ', f('normal', mu='K600_daily_mu', sigma='K600_daily_sigma')
        ) else if(pool_K600 %in% c('linear','binned')) s(
          'K600_daily', M('j'), ' ~ ', f('normal', mu=paste0('K600_daily_pred', M('j')), sigma='K600_daily_sigma')
        )
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
          condloop(
            iter='k', over=switch(pool_K600, linear='1:2', binned='1:b'),
            s('K600_daily_beta', M('k'), ' ~ ', f('normal', mu=paste0('K600_daily_beta_mu', M('k')), sigma=paste0('K600_daily_beta_sigma', M('k'))))
          )
        ),
        s('K600_daily_sigma ~ ', f('gamma', shape='K600_daily_sigma_shape', rate='K600_daily_sigma_rate'))
      )
    ),
    
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
    ode_method=c('pairmeans','Euler'),
    deficit_src=c('DO_mod','DO_obs'),
    engine=c('stan','jags'),
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
