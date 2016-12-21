#' Extract the daily parameter names from a metabolism model.
#' 
#' A function in the metab_model_interface. Returns vectors of the required and 
#' optional daily metabolism parameters for the model.
#' 
#' @inheritParams get_param_names
#' @param metab_model For get_param_names, metab_model may also be a character
#'   model name such as those created by \code{\link{mm_name}}.
#' @export
#' @examples
#' get_param_names(mm_name('mle', GPP_fun='satlight'))
#' get_param_names(mm_name('bayes'))
#' get_param_names(mm_name('Kmodel'))
#' get_param_names(mm_name('night'))
#' get_param_names(mm_name('sim'))
#' @family metab_model_interface
#' @family get_param_names
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

#' Extract the daily parameter names from a metabolism model.
#' 
#' A function in the metab_model_interface. Returns vectors of the required and 
#' optional daily metabolism parameters for the model.
#' 
#' @inheritParams get_param_names
#' @export
#' @examples 
#' dat <- data_metab()
#' get_param_names(metab(mm_name('mle', ER_fun='q10temp'), data=dat))
#' get_param_names(metab(mm_name('mle'), data=dat))
#' get_param_names(metab(mm_name('mle'), data=dat))
#' get_param_names(metab(mm_name('mle'), data=dat))
#' get_param_names(metab(mm_name('mle'), data=dat))
#' @family metab_model_interface
#' @family get_param_names
get_param_names.metab_model <- function(metab_model) {
  get_param_names(get_specs(metab_model)$model_name)
}
