#' @include metab_model-class.R
NULL

#' Collect the daily fitted parameters needed to predict GPP, ER, D, and DO
#' 
#' Returns a data.frame of parameters needed to predict GPP, ER, D, and DO
#' 
#' @inheritParams get_params
#' @return A data.frame of fitted parameters, as for the generic 
#'   \code{\link{get_params}}.
#' @export
#' @family get_params
get_params.metab_model <- function(
  metab_model, date_start=NA, date_end=NA, 
  uncertainty=c('sd','ci','none'), messages=TRUE, fixed=c('none','columns','stars'), 
  ..., attach.units=FALSE) {
  
  # process arguments
  uncertainty <- match.arg(uncertainty)
  fixed <- match.arg(fixed)
  
  # return quickly if the model is just a shell
  sp <- get_specs(metab_model)
  if(!exists('model_name', sp)) return(NULL)
  
  # build the dDOdt function in order to pull out the metab.needs
  features <- mm_parse_name(get_specs(metab_model)$model_name)
  if(features$type == 'Kmodel') {
    metab.needs <- c('K600.daily')
    metab.optional <- c()
  } else {
    egdat <- get_data(metab_model) %>% v() %>% .[which(complete.cases(.))[1:2], ]
    dDOdt <- create_calc_dDOdt(
      egdat, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
      ER_fun=features$ER_fun, deficit_src=features$deficit_src)
    metab.needs <- environment(dDOdt)$metab.needs
    metab.optional <- c('DO.mod.1') # maybe should embed this in create_calc_DO?
  }
  metab.all <- union(metab.needs, metab.optional)
  . <- '.dplyr.var'
  metab.search <- c(paste0(c('date','warnings','errors'),'$'), metab.all) %>%
    paste0('^', .) %>%
    paste0(collapse='|')
  
  # extract the daily parameters plus whatever else is daily (sds, gradients, etc.) 
  fit <- get_fit(metab_model)
  ddat <- get_data_daily(metab_model)
  
  # make sure we've got everything we need
  if(length(missing.metabs <- metab.needs[!metab.needs %in% union(names(fit), names(ddat))]) > 0) {
    stop(paste0("can't find metabolism parameter", if(length(missing.metabs)>1) "s", " ", paste0(missing.metabs, collapse=', ')))
  }
  
  # combine all daily values into one data.frame
  if(!is.null(fit) && !is.null(ddat) && nrow(ddat) > 0) {
    pars <- full_join(fit, ddat, by='date', copy=TRUE) 
  } else {
    if(!is.null(fit)) 
      pars <- fit
    else if(!is.null(ddat))
      pars <- ddat
    else
      return(NULL) # nothing available
  }
  pars <- pars %>%
    mm_filter_dates(date_start=date_start, date_end=date_end) %>% 
    { .[grep(metab.search, names(.), value=TRUE)] }
  
  # track provenance of each metab parameter. if any variables were available in
  # both x and y forms, combine them to minimize NAs
  metab.fit <- names(fit) %>% {.[. %in% metab.all]}
  metab.ddat <- names(ddat) %>% {.[. %in% metab.all]}
  metab.both <- intersect(metab.fit, metab.ddat)
  metab.either <- union(metab.fit, metab.ddat)
  for(a in metab.either) {
    if(a %in% metab.both) {
      a.x <- paste0(a,'.x')
      a.y <- paste0(a,'.y')
      pars[[a]] <- coalesce(pars[[a.x]], pars[[a.y]])
      pars[[paste0(a,'.fixed')]] <- coalesce(ifelse(is.na(pars[[a.x]]), NA, FALSE), ifelse(is.na(pars[[a.y]]), NA, TRUE))
    } else {
      pars[[paste0(a,'.fixed')]] <- a %in% metab.ddat
    }
  }
  
  # identify what we actually have, in the order we want it
  metab.out <- metab.all[metab.all %in% names(pars)]
  
  # add uncertainty columns if requested
  if(uncertainty != 'none') {
    metab.vars <- metab.out
    metab.uncert <- matrix(paste0(rep(metab.out, each=3), rep(c('.sd','.lower','.upper'), times=length(metab.out))), nrow=3, byrow=FALSE)
    metab.out <- c(rbind(metab.out, metab.uncert)) %>% { .[. %in% names(pars)]}
  }
  
  # add .fixed columns to the list of exported columns if requested
  if(fixed %in% c('columns','stars')) {
    for(a in metab.either) {
      add.after <- tail(grep(paste0('^', a), metab.out), 1)
      metab.out <- append(metab.out, paste0(a,'.fixed'), after=add.after)
    }
  }
  
  # select and order those columns of pars that match metab.needs, 
  # metab.optional, or other columns we've added. useful to order now because 
  # mm_sd_to_ci will swap columns in place
  params <- pars[c('date', metab.out)]
  
  # convert sds to CIs if requested
  if(uncertainty == 'sd') {
    extra.cols <- grep('\\.lower$|\\.upper$', names(params))
    if(length(extra.cols) > 0) params <- params[-extra.cols]
  } else if(uncertainty == 'ci') {
    # use existing .lower and .upper cols if available
    for(mv in metab.vars) {
      if(all(paste0(mv,c('.lower','.upper')) %in% names(params))) {
        extra.cols <- grep(paste0(mv,'\\.sd$'), names(params))
        if(length(extra.cols) > 0) params <- params[-extra.cols]
      }
    }
    # convert any remaining .sd cols to .lower and .upper parametrically
    params <- mm_sd_to_ci(params)
  }
  
  # convert .fixed columns to stars if requested (do this after mm_sd_to_ci b/c converts to character)
  if(fixed == 'stars') {
    params <- bind_cols(select(params, date), format.data.frame(select(params, -date)))
    for(a in metab.either) {
      params[[a]] <- paste0(params[[a]], ifelse(params[[paste0(a,'.fixed')]], '*', ' '))
      params[[paste0(a,'.fixed')]] <- NULL
    }
  }
  
  # attach warnings and errors if requested
  if(messages && exists('date', pars) && any(exists(c('warnings','errors'), pars))) {
    messages <- pars[c('date','warnings','errors') %>% { .[. %in% names(pars)] }]
    pretty_print_ddat
    params <- left_join(params, messages, by='date', copy=TRUE)
  }
  
  # attach units if requested and available in mm_data
  if(attach.units) {
    param.units <- get_units(mm_data())[names(params)]
    params <- u(params, param.units)
  }
  
  # return
  params
}
