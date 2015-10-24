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
  warning(mtb_funs)
  
  if(!all(mtb_funs %in% c("metab_bayes","metab_Kvpred", "metab_mle", "metab_night", "metab_sim")) || length(mtb_funs) != 5)
    stop("metab_funs() output is out of date")
  
  data.frame(
    metab_fun=c("metab_bayes","metab_Kvpred", "metab_mle", "metab_night", "metab_sim"),
    method=c("Bayesian MCMC","regression of daily values","maximum likelihood estimation","nighttime regression","data simulation")
  )
  
}