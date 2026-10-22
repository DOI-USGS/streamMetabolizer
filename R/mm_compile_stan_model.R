#' Compile a Stan model, or load it from a cached .stanrds file
#'
#' Shared by \code{runstan_bayes} (one-station) and \code{metab_bayes_2s}
#' (two-station): compiles \code{model_path} via \code{rstan::stan_model()}
#' and caches the result alongside it as \code{<model>.stanrds}, so a
#' second call against the same file loads the cached object instead of
#' recompiling. The cache is invalidated whenever \code{model_path} is
#' newer than the cached \code{.stanrds} file.
#'
#' @param model_path the Stan model file to compile, as a full file path
#' @param verbose logical. give status messages?
#' @return a list with \code{stan_mobj} (the compiled Stan model object),
#'   \code{compile_time} (a \code{proc_time}, zero if loaded from cache),
#'   and \code{compile_log} (character vector of compiler output, or
#'   \code{NULL} if loaded from cache)
#' @keywords internal
mm_compile_stan_model <- function(model_path, verbose=FALSE) {

  # stan_model()/sampling() can't find their own function
  # cpp_object_initializer() unless the namespace is loaded. requireNamespace
  # is somehow not doing this. Thoughts (not solution):
  # https://stat.ethz.ch/pipermail/r-devel/2014-September/069803.html
  if(!suppressPackageStartupMessages(require(rstan))) {
    stop("the rstan package is required for Stan MCMC models")
  }

  # use auto_write=TRUE to recompile if needed, or load from existing .rds file
  # without recompiling if possible
  compile_time <- system.time({})
  compile_log <- NULL
  mobj_path <- gsub('.stan$', '.stanrds', model_path)
  if(!file.exists(mobj_path) || file.info(mobj_path)$mtime < file.info(model_path)$mtime) {
    if(verbose) message("compiling Stan model")
    compile_time <- system.time({
      compile_log <- capture.output({
        stan_mobj <- rstan::stan_model(file=model_path, auto_write=TRUE)
      }, type=c('output'), split=verbose)
    })
    rm(stan_mobj)
    gc() # this humble line saves us from many horrible R crashes
    autowrite_path <- gsub('.stan$', '.rds', model_path)
    if(!file.exists(autowrite_path)) autowrite_path <- gsub('.stan$', '.rda', model_path) # for backwards compatibility with rstan < 2.13
    if(!file.exists(autowrite_path)) autowrite_path <- file.path(tempdir(), basename(autowrite_path))
    if(!file.exists(autowrite_path)) {
      warning('could not find saved rds model file')
    } else {
      tryCatch({
        file.copy(autowrite_path, mobj_path, overwrite=TRUE)
        file.remove(autowrite_path)
      }, error=function(e) {
        warning('could not copy Stan rds to .stanrds file: ', e$message)
        mobj_path <- autowrite_path
      })
    }
  } else {
    if(verbose) message("loading pre-compiled Stan model")
  }
  stan_mobj <- readRDS(mobj_path)

  list(stan_mobj=stan_mobj, compile_time=compile_time, compile_log=compile_log)
}
