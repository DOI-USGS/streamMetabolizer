#' \code{specs_mle_obserr} - maximum likelihood estimation of GPP, ER, and 
#' optionally also K600, assuming error in [only] the observation of DO.obs
#' 
#' @rdname specs_mle
#'   
#' @inheritParams negloglik_1ply
#' @param GPP_init the inital value of daily GPP to use in the NLM fitting 
#'   process
#' @param ER_init the inital value of daily ER to use in the NLM fitting process
#' @param K600_init the inital value of daily K600 to use in the NLM fitting 
#'   process. Ignored if K600 is supplied in data_daily, except for those dates 
#'   where K600 is NA. If there are any such dates, K600_init must have a 
#'   numeric (non-NA) value, as this will be used to estimate K600 for those 
#'   dates.
#' @export
#' @family model_specs
specs_mle_obserr <- function(
  calc_DO_fun=calc_DO_mod, ODE_method="pairmeans", # inheritParams negloglik_1ply
  GPP_init=3, ER_init=-5, K600_init=10
) {
  list(
    calc_DO_fun=calc_DO_fun,
    ODE_method=ODE_method,
    GPP_init=GPP_init,
    ER_init=ER_init,
    K600_init=K600_init
  )
}