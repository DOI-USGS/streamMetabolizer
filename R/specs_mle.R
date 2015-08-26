#' Define the parameters
#' 
#' @param GPP.init the inital value of daily GPP to use in the NLM fitting
#'   process
#' @param ER.init the inital value of daily ER to use in the NLM fitting process
#' @param K600.init the inital value of daily K600 to use in the NLM fitting
#'   process
metab_spec_mle <- function(
  GPP.init=3, ER.init=-5, K600.init=5
) {
  list(
    GPP.init=GPP.init,
    ER.init=ER.init,
    K600.init=K600.init
  )
}