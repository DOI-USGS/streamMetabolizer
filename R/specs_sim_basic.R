#' Simulation of DO.obs from input data and metabolism
#' 
#' @inheritParams specs_all
#'   
#' @export
specs_sim_basic <- function(
 
  err.obs.sigma=0.1, 
  err.obs.phi=0, 
  err.proc.sigma=0, 
  err.proc.phi=0,
  ODE_method="pairmeans",
  sim.seed=NA
  
) {
  
  as.list(environment())
  
}