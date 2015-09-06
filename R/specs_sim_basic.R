#' \code{specs_sim_basic} - model for simulating DO.obs from input data and
#' metabolism
#' 
#' @rdname specs_sim
#'
#' @export
#' @family model_specs
specs_sim_basic <- function(
  obs.err.sd, 
  obs.err.phi, 
  proc.err.sd, 
  proc.err.phi
) {
  list(
    obs.err.sd=obs.err.sd, 
    obs.err.phi=obs.err.phi, 
    proc.err.sd=proc.err.sd, 
    proc.err.phi=proc.err.phi
  )
}