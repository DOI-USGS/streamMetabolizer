#' List the metabolism model fitting functions
#' 
#' @import dplyr
#' @export
metab_funs <- function() {
  # list the metab_ functions
  . <- '.dplyr.var'
  mtb_funs <- ls("package:streamMetabolizer") %>%
    grep("^metab_", ., value=TRUE) %>% 
    grep("model*|funs", ., invert=TRUE, value=TRUE)
  
  known_funs <- c("metab_bayes","metab_Kmodel", "metab_mle", "metab_night", "metab_sim")
  if(!all(mtb_funs %in% known_funs) || length(mtb_funs) != 5)
    stop("metab_funs() output is out of date")
  
  data.frame(
    metab_fun=known_funs,
    method=c("Bayesian MCMC","statistical model of daily K values","maximum likelihood estimation","nighttime regression","data simulation")
  )
  
}