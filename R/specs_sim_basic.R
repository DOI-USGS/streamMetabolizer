#' \code{specs_sim_basic} - model for simulating DO.obs from input data and 
#' metabolism
#' 
#' @rdname specs_sim
#'   
#' @inheritParams calc_DO_mod_w_sim_error
#' @param sim.seed NA to specify that each call to predict_DO should generate 
#'   new values, or an integer, as in the \code{seed} argument to 
#'   \code{\link{set.seed}}, specifying the seed to set before every execution
#'   of predict_DO
#'   
#' @export
#' @family model_specs
specs_sim_basic <- function(
  err.obs.sigma=0.1, 
  err.obs.phi=0, 
  err.proc.sigma=0, 
  err.proc.phi=0,
  sim.seed=NA
) {
  list(
    err.obs.sigma=err.obs.sigma, 
    err.obs.phi=err.obs.phi, 
    err.proc.sigma=err.proc.sigma, 
    err.proc.phi=err.proc.phi,
    sim.seed=sim.seed
  )
}