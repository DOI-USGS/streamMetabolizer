#' Identify a model_specs function by model traits
#' 
#' @param metab_fun The model fitting function, which defines the method 
#'   (metab_mle = Maximum Likelihood Estimation, metab_bayes = Bayesian MCMC, 
#'   etc.)
#' @param software The modeling software used, if any, in addition to R (mainly
#'   applies to Bayesian methods: JAGS, Stan, etc.)
#' @param pool_days Does the model use pooling to combine information among 
#'   days?
#' @param err_obs_iid Does the model accommodate IID (independent and 
#'   identically distributed) observation error?
#' @param err_proc_iid Does the model accommodate IID process error? Process 
#'   error causes autocorrelation in total prediction error but needn't be 
#'   autocorrelated itself. If autocorrelation parameters are not defined for 
#'   the process errors, we consider the process errors to be IID.
#' @param err_proc_acor Does the model accommodate autocorrelated process error?
#'   If the model includes autocorrelation parameters for process error, even if
#'   the autocorrelation coefficient might ultimately be given a value of 0, we 
#'   consider the model to accommodate autocorrelated error.
#' @import dplyr
#' @export
#' @examples
#' specs_funs()
#' specs_funs(metab_fun="metab_mle", err_obs_iid=TRUE, err_proc_iid=FALSE)$specs_fun
#' dplyr::filter(specs_funs(), err_proc_iid | err_proc_acor)
specs_funs <- function(metab_fun=c("metab_bayes","metab_mle","metab_night","metab_sim"), 
                       software=c("jags","stan",NA), pool_days=c(TRUE,FALSE),
                        err_obs_iid=c(TRUE,FALSE), err_proc_iid=c(TRUE,FALSE), err_proc_acor=c(TRUE,FALSE)) {
  models <- data.frame(
    specs_fun=paste0(
      "specs_",
      c("bayes_jags_nopool_obserr",
        "bayes_jags_nopool_procobserr",
        "bayes_stan_nopool_obserr",
        "bayes_stan_nopool_procacoriiderr",
        "bayes_stan_nopool_procobserr",
        "mle_obserr",
        "mle_procerr",
        "night_basic",
        "sim_basic"
      )),
    metab_fun=paste0("metab_", c(rep("bayes",5), rep("mle",2), "night", "sim")),
    software=c(rep("jags",2), rep("stan",3), rep(NA,4)), 
    pool_days=    rep(F, 9), 
    err_obs_iid=  c(T,T,T,F,T,T,F,T,T), 
    err_proc_iid= c(F,F,F,T,F,F,T,F,T),
    err_proc_acor=c(F,T,F,T,T,F,F,F,T), 
    stringsAsFactors = FALSE)
  
  metab_fun_ <- software_ <- pool_days_ <- err_obs_iid_ <- err_proc_iid_ <- err_proc_acor_ <- '.dplyr.var'
  models %>% 
    setNames(paste0(names(models), "_")) %>%
    dplyr::filter(
      tolower(metab_fun_) %in% tolower(metab_fun),
      tolower(software_) %in% tolower(software),
      pool_days_ %in% pool_days,
      err_obs_iid_ %in% err_obs_iid,
      err_proc_iid_ %in% err_proc_iid,
      err_proc_acor_ %in% err_proc_acor) %>%
    setNames(names(models))
  
}