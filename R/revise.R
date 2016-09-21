#' Change or add named elements of a list
#' 
#' Primary use case is for revising a list of specifications as originally 
#' created by specs()
#' 
#' @param specs a list of specifications to revise
#' @param ... named values to replace in or add to \code{specs}
#' @param delete a vector of names of elements to remove from the specs list
#' @export
#' @examples
#' sp <- specs(mm_name('bayes'))
#' sp <- revise(sp, 
#'   model_name='b_np_oipi_tr_plrckm_mynewmodel.stan',
#'   params_in=c(params_in,'my_new_param'), my_new_param=4)
revise <- function(specs, ..., delete) {
  # evaluate ... in the context of specs
  args <- {
    attach(specs, warn.conflicts=FALSE, name='revisespecstemp')
    on.exit(detach('revisespecstemp'))
    list(...)
  }
  
  # require names
  if(is.null(names(args)) || any(names(args) == '')) {
    stop("all arguments in ... must be named")
  }
  # add/change specs
  for(a in names(args)) {
    specs[[a]] <- args[[a]]
  }
  # delete specs
  if(!missing(delete)) {
    specs[delete] <- NULL
  }
  # return
  specs
}
