#' MLE estimation with process error only
#' 
#' Maximum likelihood estimation of GPP, ER, and optionally also K600, assuming 
#' error in [only] the process of reaching DO.obs[t+1] from DO.obs[t]
#' 
#' @inheritParams specs_all
#'   
#' @export
specs_mle_procerr <- function(
  
  # inheritParams negloglik_1ply
  calc_DO_fun=calc_DO_mod_by_diff, 
  ODE_method="pairmeans",
  
  # inheritParams specs_mle_obserr
  GPP_init=3, 
  ER_init=-5, 
  K600_init=10
  
) {
  
  as.list(environment())
  
}