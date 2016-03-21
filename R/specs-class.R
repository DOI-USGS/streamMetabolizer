#### initialize ####

setOldClass('specs')

#' Add the metabolism model specifications class to a list
#' 
#' @keywords internal
add_specs_class <- function(specs_list) {
  class(specs_list) <- c("specs", class(specs_list))
  specs_list
}

#### display ####

#' Display the specs object
#' 
#' Print a specs object to the console.
#' 
#' @param object specs list to be displayed.
setMethod(
  "show", "specs", 
  function(object) {
    print_specs(object)
  }
)

#' Display the specs object
#' 
#' Print a specs object to the console.
#' 
#' @param x specs list to be displayed.
#' @param ... additional arguments passed to inner functions
#' @export
print.specs <- function(x, ...) {
  print_specs(x, ...)
}

#' Display the specs object
#' 
#' Print a specs object to the console.
#' 
#' @param object specs list to be displayed.
#' @param header line to be catted at start of printout
#' @param prefix text to prepend to the start of each line that follows the
#'   header
#' @keywords internal
print_specs <- function(object, header="Model specifications:\n", prefix="  ") {
  cat(header)
  for(spec in names(object)) {
    spec_char <- tryCatch(
      if(is.null(object[[spec]])) {
        'NULL' 
      } else {
        paste0(as.character(object[[spec]]), collapse=", ")
      }, 
      error=function(e) paste0(paste0(class(object[[spec]]), collapse=","),"; see element [['",spec,"']] for details"))
    if(nchar(spec_char) > 100) spec_char <- paste0(substr(spec_char, 1, 100), "...")
    cat(paste0(prefix, spec, "\t", spec_char, "\n"))
  }
}