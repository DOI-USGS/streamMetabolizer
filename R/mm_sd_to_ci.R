#' Convert SD columns into CI columns in a data.frame
#' 
#' Convert data with var and var.sd columns into data with var, var.lower, and 
#' var.upper columns
#' 
#' @param data a data.frame with 1+ pairs of columns named var and var.sd 
#'   (where var can be anything)
#' @param alpha the desired significance level described by the confidence 
#'   interval
#' @import dplyr
#' @keywords internal
mm_sd_to_ci <- function(data, alpha=0.05) {
  
  # identify the pairs of columns
  . <- '.dplyr.var'
  sd.cols <- grep('\\.sd$', names(data), value=TRUE)
  var.cols <- substring(sd.cols, 1, nchar(sd.cols)-3) %>% { .[. %in% names(data)] } # could be fewer than sd.cols
  if(length(sd.cols) > 0) sd.cols <- paste0(var.cols, '.sd') # only include those cols that both end in '.sd' and have a paired var.cols value
  sd.pos <- match(sd.cols, names(data))
  
  # replace each var.sd col with var.lower and var.upper cols
  crit <- -qnorm(alpha/2)
  dat.list <- as.list(data)
  for(i in rev(seq_along(var.cols))) {
    est <- dat.list[[var.cols[i]]]
    sd <- dat.list[[sd.cols[i]]]
    dat.list <- append(
      dat.list, 
      setNames(list(
        lower = est - crit * sd,
        upper = est + crit * sd), paste0(var.cols[i], c('.lower','.upper'))),
      after=sd.pos[i])
    dat.list[sd.pos[i]] <- NULL
  }
  
  as.data.frame(dat.list, stringsAsFactors=FALSE)
}
