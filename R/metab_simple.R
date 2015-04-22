#'@title Basic metabolism model fitting function
#'@description Fits a model to estimate GPP and ER from input data on DO, 
#'  temperature, light, etc.
#'  
#'@param DO.obs Vector of dissolved oxygen concentration observations, \eqn{mg 
#'  O[2] L^{-1}}{mg O2 / L}
#'@param DO.sat Vector of dissolved oxygen saturation values based on water 
#'  temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models 
#'  (e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param PAR Vector of photosynthetically active radiation in \eqn{\mu mol\ 
#'  m^{-2} s^{-1}}{micro mols / m^2 / s}
#'@param temp.water Vector of water temperatures in \eqn{^{\circ}C}{degrees C}. Used in
#'  scaling respiration with temperature
#'@param ... additional arguments
#'@return A data.frame with columns corresponding to components of metabolism 
#'  \describe{
#'  
#'  \item{GPP}{numeric estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} 
#'  d^{-1}}{mg O2 / L / d}}
#'  
#'  \item{R}{numeric estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 /
#'  L / d}}
#'  
#'  \item{NEP}{numeric estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} 
#'  d^{-1}}{mg O2 / L / d}}
#'  
#'  }
#'  
#'  Use \link{attributes} to access more model output:
#'  
#'  \item{params}{parameter estimates}
#'  
#'@details The model has inputs and parameters
#'  
#'@references 
#'@seealso 
#'@author Alison Appling, Jordan Read; modeled on lakeMetabolizer
#'@examples
#'  metab_simple(date.time=rep(as.Date("2015-04-15"), 24), DO.obs=1:24, DO.sat=rep(14, 24), PAR=sin((1:24)*pi/24)^8, temp.water=15)
#'@export
metab_simple <- function(date.time, DO.obs, DO.sat, k.gas, PAR, temp.water, ...) {
  set.seed(3000)
  GPP <- abs(rnorm(length(unique(as.Date(date.time))), 4, 1))
  ER <- -1.6*GPP
  data.frame(GPP=GPP, ER=ER, NEP=GPP-ER)
}