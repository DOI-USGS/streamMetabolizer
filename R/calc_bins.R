#' Assign continuous values in a vector to discrete bins
#' 
#' Assigns each value in \code{vec} a new, discrete value corresponding to a 
#' bin. This function provides one interface to the functions 
#' `base::cut`, `ggplot2::cut_interval`, and `ggplot2::cut_number`.
#' 
#' @param vec the numeric vector whose values should be binned. 
#'   log(discharge.daily) is a good candidate when using this function for 
#'   pooling of K600 values.
#' @param method a single character string indicating the automated bin
#'   selection method to use
#' @param bounds if method=='bounds', a numeric vector of bin boundaries
#' @param \dots other arguments (e.g. \code{n}, \code{width}) passed to the 
#'   ggplot function corresponding to the value of cuts, if cuts is a character 
#'   (otherwise ignored)
#' @importFrom unitted v
#' @export
#' @examples
#' ln.disch <- log(rlnorm(100))
#' 
#' # for use in setting specs
#' brks <- calc_bins(ln.disch, 'width', width=0.8)$bounds
#' specs('b_Kb_oipi_tr_plrckm.stan', K600_lnQ_nodes_centers=brks)
#' 
#' # variations
#' 
#' # by 'number' method
#' bins_num <- calc_bins(ln.disch, 'number', n=5)
#' df_num <- data.frame(t=1:length(ln.disch), vec=ln.disch, bin=bins_num$names[bins_num$vec])
#' table(bins_num$vec)
#' 
#' # by 'interval' method
#' bins_int <- calc_bins(ln.disch, 'interval', n=5)
#' df_int <- data.frame(t=1:length(ln.disch), vec=ln.disch, bin=bins_int$names[bins_int$vec])
#' table(bins_int$vec)
#' 
#' # by 'width' method
#' bins_wid <- calc_bins(ln.disch, 'width', width=0.2, boundary=0)
#' df_wid <- data.frame(t=1:length(ln.disch), vec=ln.disch, bin=bins_wid$names[bins_wid$vec])
#' table(bins_wid$vec)
#' 
#' # choose your own arbitrary breaks
#' bins_arb <- calc_bins(ln.disch, bounds=seq(-4,4,by=1))
#' df_arb <- data.frame(t=1:length(ln.disch), vec=ln.disch, bin=bins_arb$names[bins_arb$vec])
#' table(bins_arb$vec)
#' \dontrun{
#' library(ggplot2)
#' ggplot(df_num, aes(x=t, y=vec, color=bin)) + geom_point() + 
#'   geom_hline(data=as.data.frame(bins_num['bounds']), aes(yintercept=bounds))
#' ggplot(df_int, aes(x=t, y=vec, color=bin)) + geom_point() +
#'   geom_hline(data=as.data.frame(bins_int['bounds']), aes(yintercept=bounds))
#' ggplot(df_wid, aes(x=t, y=vec, color=bin)) + geom_point() +
#'   geom_hline(data=as.data.frame(bins_wid['bounds']), aes(yintercept=bounds))
#' ggplot(df_arb, aes(x=t, y=vec, color=bin)) + geom_point() +
#'   geom_hline(data=as.data.frame(bins_arb['bounds']), aes(yintercept=bounds))
#' }
calc_bins <- function(vec, method=c('bounds','interval','number','width'), ..., bounds) {
  
  method <- match.arg(method)
  if(method != 'bounds') {
    
    if(!requireNamespace('ggplot2', quietly=TRUE)) {
      stop("need ggplot2 to calculate discharge bins when is.character(method). ",
           "either install ggplot2 or switch to a numeric vector for method")
    }
    # run once with high dig.lab to parse the breaks from levels(cutvals) as numeric
    cutvals <- switch(
      method,
      interval = ggplot2::cut_interval(v(vec), ..., dig.lab=20),
      number = ggplot2::cut_number(v(vec), ..., dig.lab=20),
      width = ggplot2::cut_width(v(vec), ...)) # dig.lab is unavailable for width
    bounds <- levels(cutvals) %>%
      strsplit('\\[|\\(|\\]|,') %>%
      lapply(function(lev) as.numeric(lev[2:3])) %>%
      unlist() %>%
      unique()
    bounds[c(1,length(bounds))] <- bounds[c(1,length(bounds))] + c(-1e-10, 1e-10) # add a fudge factor to the outer bounds
    
  } # else we already know bounds precisely because we're supplying them
  
  # final run with precise numeric breaks that were chosen
  binned <- cut(vec, breaks=bounds, ordered_result=TRUE)
  
  list(vec=as.numeric(binned), bounds=bounds, names=levels(binned))
}
