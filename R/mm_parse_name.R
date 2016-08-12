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
#' @import dplyr
#' @importFrom stats na.omit
#' @examples
#' mm_parse_name(c(mm_name('mle'), mm_name('night'), mm_name('bayes')))
#' mm_parse_name(c(mm_name('mle'), mm_name('night'), mm_name('bayes')), keep_name=TRUE)
#' @export
mm_parse_name <- function(model_name, keep_name=FALSE) {

  # define function that gets used to parse prk_terms
  match_or_NA <- function(key, pairs) { 
    matches <- c(unname(na.omit(key[pairs]))) 
    if(length(matches) == 0) {
      'NA'
    } else if(length(matches) > 1) {
      stop('found too many matches in PRK terms') 
    } else { matches }
  }

  # parse the name
  parsed <- strsplit(basename(model_name), "_|\\.")
  sapply(1:length(parsed), function(pnum) if(length(parsed[[pnum]]) <= 5) stop('missing one or more pieces in name: ', model_name[pnum]))
  type <- unname(c(b='bayes', m='mle', n='night', K='Kmodel', s='sim')[sapply(parsed, `[`, 1)])
  pool_K600 <- unname(c(np='none', Kn='normal', Kl='linear', Kb='binned', Kc='complete')[sapply(parsed, `[`, 2)])
  err_obs_iid <- grepl('oi', sapply(parsed, `[`, 3))
  err_proc_acor <- grepl('pc', sapply(parsed, `[`, 3))
  err_proc_iid <-  grepl('pi', sapply(parsed, `[`, 3))
  ode_method <- unname(
    c(Eu='Euler', pm='pairmeans', tr='trapezoid', r2='rk2', o1='lsoda', o2='lsode', o3='lsodes', 
      o4='lsodar', o5='vode', o6='daspk', o7='euler', eu='euler', o8='rk4', o9='ode23', o10='ode45', o11='radau', 
      o12='bdf', o13='bdf_d', o14='adams', o15='impAdams', o16='impAdams_d')[sapply(parsed, `[`, 4)])
  prk_terms <- bind_rows(lapply(parsed, function(parsed1) {
    prk_term <- parsed1[5]
    prk_pairs <- if(nchar(prk_term)==0) c() else sapply(seq(2, nchar(prk_term), by=2), function(pos) substring(prk_term, pos-1, pos))
    data.frame(
      deficit_src=match_or_NA(c(km='DO_mod', ko='DO_obs', kf='DO_obs_filter'), prk_pairs),
      ER_fun=match_or_NA(c(rc='constant', rq='q10temp'), prk_pairs),
      GPP_fun=match_or_NA(c(pl='linlight', ps='satlight'), prk_pairs),
      stringsAsFactors=FALSE)
  }))
  GPP_fun <- prk_terms$GPP_fun
  ER_fun <- prk_terms$ER_fun
  deficit_src <- prk_terms$deficit_src
  engine <- sapply(parsed, function(vec) vec[length(vec)]) # the last one - leaves room for custom name endings before the suffix
  
  # combine the parsed pieces into a data.frame
  df <- data.frame(
    model_name=model_name,
    type=type,
    pool_K600=pool_K600,
    err_obs_iid=err_obs_iid,
    err_proc_acor=err_proc_acor,
    err_proc_iid=err_proc_iid,
    ode_method=ifelse(is.na(ode_method), 'NA', ode_method),
    GPP_fun=ifelse(is.na(GPP_fun), 'NA', GPP_fun),
    ER_fun=ifelse(is.na(ER_fun), 'NA', ER_fun),
    deficit_src=ifelse(is.na(deficit_src), 'NA', deficit_src),
    engine=ifelse(is.na(engine), 'NA', engine), 
    stringsAsFactors=FALSE)
  
  if(!keep_name) df$model_name <- NULL
  
  df
}