#' MLE estimation with observation error only
#' 
#' Maximum likelihood estimation of GPP, ER, and optionally also K600, assuming 
#' error in [only] the observation of DO.obs
#' 
#' @inheritParams specs_all
#'   
#' @export
specs_mle_obserr <- function(
  
  # inheritParams negloglik_1ply
  calc_DO_fun=calc_DO_mod, 
  ODE_method="pairmeans",
  
  GPP_init=3, 
  ER_init=-5, 
  K600_init=10
  
) {
  
  as.list(environment())
  
}