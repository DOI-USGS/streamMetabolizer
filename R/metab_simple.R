#' @title Basic metabolism model fitting function
#' @description Fits a model to estimate GPP and ER from input data on DO, 
#'   temperature, light, etc.
#'   
#' @param data data.frame with columns \itemize{ \item{ DO.obs Vector of 
#'   dissolved oxygen concentration observations, \eqn{mg O[2] L^{-1}}{mg O2 / 
#'   L}} \item{ DO.sat Vector of dissolved oxygen saturation values based on 
#'   water temperature. Calculate using \link{o2.at.sat}} \item{ k.gas Vector of
#'   kGAS values calculated from any of the gas flux models (e.g., 
#'   \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}} \item{ PAR 
#'   Vector of photosynthetically active radiation in \eqn{\mu mol\ m^{-2} 
#'   s^{-1}}{micro mols / m^2 / s}} \item{ temp.water Vector of water 
#'   temperatures in \eqn{^{\circ}C}{degrees C}. Used in scaling respiration 
#'   with temperature} }
#' @param ... additional arguments
#' @return A data.frame with columns corresponding to components of metabolism 
#'   \describe{ \item{GPP}{numeric estimate of Gross Primary Production, \eqn{mg
#'   O_2 L^{-1} d^{-1}}{mg O2 / L / d}} \item{R}{numeric estimate of
#'   Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}} \item{NEP}{numeric
#'   estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L
#'   / d}} }
#'   
#' @details The model has inputs and parameters
#'   
#' @author Alison Appling, Jordan Read; modeled on lakeMetabolizer
#' @examples
#'  metab_simple(data=data.frame(
#'    date.time=rep(as.Date("2015-04-15"), 24), DO=1:24, 
#'    DO.deficit=rep(2, 24), depth=sin((1:24)*pi/24)^8, k.DO=15))
#'  \dontrun{
#'  metab_simple(data=data.frame(empty="shouldbreak"))
#'  }
#' @export
metab_simple <- function(data, ...) {
  
  # Check data for correct column names
  expected_colnames <- c("date.time","DO","DO.deficit","depth","k.DO")
  if(!all(expected_colnames %in% names(data))) {
    stop(paste0("data must contain (at least) columns with the names ", paste0(expected_colnames, collapse=", ")))
  }
  
  #' Return the likelihood value for a given set of parameters and observations
  #' 
  #' From ?nlm, this function should be "the function to be minimized, returning
  #' a single numeric value. This should be a function with first argument a 
  #' vector of the length of p followed by any other arguments specified by the 
  #' ... argument."
  #' 
  #' @param params a vector of length 2, where the first element is GPP and the second element is ER (both mg/L/d)
  #' @param 
#   onestation_likelihood <- function(params, DO, DO.deficit, oxy, light, z, bp, ts, K) {
#     
#     GPP <- params[1]
#     ER <- params[2]
#     
#     DO.modeled <- numeric(nrow(data))
#     DO.modeled[1] <- oxy[1]
#     for (i in 2:length(DO)) {
#       DO.modeled[i] <- 
#         DO.modeled[i-1] +
#         (GPP/z)*(light[i]/sum(light)) + 
#         ER * ts / z + 
#         Kcor(temp[i],K)*ts*DO.deficit 
#     }
#   
#     ##below is MLE calculation
#     sqdiff<-(oxy-metab)^2 
#     length(oxy)*(log(((sum(sqdiff)/length(oxy))^0.5)) +0.5*log(6.28))   + ((2*sum(sqdiff)/length(oxy))^-1)*sum(sqdiff)
#   }
  
  
  # Calculate metabolism by non linear minimization of an MLE function
#   river.mle <- nlm(mlefunction, p=c(3,-5), oxy=oxy, z=z,temp=temp,light=light, bp=bp, ts=ts, K=K)
  
#   b<-c(river.mle$estimate[1], river.mle$estimate[2],river.mle$estimate[3],  river.mle$minimum[1])
  
  metab_model()
}