#' \code{specs_mle_procerr} - maximum likelihood estimation of GPP, ER, and 
#' optionally also K600, assuming error in [only] the process of reaching
#' DO.obs[t+1] from DO.obs[t]
#' 
#' @rdname specs_mle
#'   
#' @inheritParams negloglik_1ply
#' @inheritParams specs_mle_obserr
#' @export
#' @family model_specs
specs_mle_procerr <- function(
  calc_DO_fun=calc_DO_mod_by_diff, ODE_method="pairmeans", # inheritParams negloglik_1ply
  GPP_init=3, ER_init=-5, K600_init=10 # inheritParams specs_mle_obserr
) {
  list(
    calc_DO_fun=calc_DO_fun,
    ODE_method=ODE_method,
    GPP_init=GPP_init,
    ER_init=ER_init,
    K600_init=K600_init
  )
}