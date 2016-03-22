#' Parse a model name into its features
#' 
#' Returns a data.frame with one column per model structure detail and one row 
#' per `model_name` supplied to this function. See \code{?\link{mm_name}} for a 
#' description of each of the data.frame columns that is returned.
#' 
#' Custom model files (for MCMC) may have additional characters after an 
#' underscore at the end of the name and before the prefix. For example, 
#' 'b_np_pcpi_eu_ko.stan' and 'b_np_pcpi_eu_ko_v2.stan' are parsed the same; the
#' _v2 is ignored by this function.
#' 
#' @seealso The converse of this function is \code{\link{mm_name}}.
#'   
#' @param model_name character: the model name
#' @param keep_name logical; should the model_name be included as a first column
#'   in the output data.frame?
#' @examples
#' mm_parse_name(c(mm_name('mle'), mm_name('night'), mm_name('bayes')))
#' mm_parse_name(c(mm_name('mle'), mm_name('night'), mm_name('bayes')), keep_name=TRUE)
#' @export
mm_parse_name <- function(model_name, keep_name=FALSE) {

  parsed <- strsplit(basename(model_name), "_|\\.")
  sapply(1:length(parsed), function(pnum) if(length(parsed[[pnum]]) <= 5) stop('missing one or more pieces in name: ', model_name[pnum]))
  type <- unname(c(b='bayes', m='mle', n='night', K='Kmodel', s='sim')[sapply(parsed, `[`, 1)])
  pool_K600 <- unname(c(np='none', Kn='normal', Kl='linear', Kb='binned')[sapply(parsed, `[`, 2)])
  err_obs_iid <- grepl('oi', sapply(parsed, `[`, 3))
  err_proc_acor <- grepl('pc', sapply(parsed, `[`, 3))
  err_proc_iid <-  grepl('pi', sapply(parsed, `[`, 3))
  ode_method <- unname(c(eu='Euler', pm='pairmeans')[sapply(parsed, `[`, 4)])
  deficit_src <- unname(c(km='DO_mod', ko='DO_obs')[sapply(parsed, `[`, 5)])
  engine <- sapply(parsed, function(vec) vec[length(vec)]) # the last one - leaves room for custom name endings before the suffix
  
  df <- data.frame(
    model_name=model_name,
    type=type,
    pool_K600=pool_K600,
    err_obs_iid=err_obs_iid,
    err_proc_acor=err_proc_acor,
    err_proc_iid=err_proc_iid,
    ode_method=ifelse(is.na(ode_method), 'NA', ode_method),
    deficit_src=ifelse(is.na(deficit_src), 'NA', deficit_src),
    engine=ifelse(is.na(engine), 'NA', engine), 
    stringsAsFactors=FALSE)
  
  if(!keep_name) df$model_name <- NULL
  
  df
}