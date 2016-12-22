#' @describeIn get_param_names This implementation is shared by many model types
#' @export
#' @examples
#' 
#' # pass in a character string:
#' get_param_names(mm_name('mle', GPP_fun='satlight'))
#' get_param_names(mm_name('bayes'))
#' get_param_names(mm_name('Kmodel'))
#' get_param_names(mm_name('night'))
#' get_param_names(mm_name('sim'))
get_param_names.character <- function(metab_model) {
  # parse the model features
  features <- mm_parse_name(metab_model)
  
  # build the dDOdt function in order to pull out the metab.needs
  if(features$type == 'Kmodel') {
    metab.needs <- c('K600.daily')
    metab.optional <- c()
  } else {
    egdat <- unitted::v(eval(formals(paste0("metab_", features$type))$data)) %>%
      bind_rows(.,.)
    dDOdt <- create_calc_dDOdt(
      egdat, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
      ER_fun=features$ER_fun, deficit_src=features$deficit_src)
    metab.needs <- environment(dDOdt)$metab.needs
    metab.optional <- c('DO.mod.1') # maybe should embed this in create_calc_DO?
  }
  list(required=metab.needs, optional=metab.optional)
}

#' @describeIn get_param_names Lets you pass in a model object rather than a
#'   character string
#' @export
#' @examples 
#' 
#' # or pass in a metab_model object:
#' dat <- data_metab()
#' get_param_names(metab(specs(mm_name('mle', ER_fun='q10temp')), data=dat))
#' get_param_names(metab(specs('bayes'), data=dat))
#' get_param_names(metab(specs('night'), data=dat))
#' get_param_names(metab(specs('Kmodel'), data=dat))
#' get_param_names(metab(specs('sim'), data=dat))
get_param_names.metab_model <- function(metab_model) {
  get_param_names(get_specs(metab_model)$model_name)
}
