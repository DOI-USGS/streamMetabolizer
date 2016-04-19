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
#' @import dplyr
#' @keywords internal
print_specs <- function(object, header="Model specifications:\n", prefix="  ") {
  # create a data.frame with a concise 1-line description of each specs element
  max_value_width <- if(length(object) > 0) {
    max(10, getOption('width') - nchar(prefix) - max(nchar(names(object))) - 1)
  } else { 
    10
  }
  specs_df <- data.frame(value=sapply(names(object), function(spec) {
    spec_char <- tryCatch(
      if(is.null(object[[spec]])) {
        'NULL' 
      } else {
        paste0(as.character(object[[spec]]), collapse=", ")
      }, 
      error=function(e) {
        paste0(paste0(class(object[[spec]]), collapse=","),"; see element [['",spec,"']] for details")
      })
    if(nchar(spec_char) > max_value_width) {
      spec_char <- paste0(substr(spec_char, 1, max_value_width-3), "...")
    }
    spec_char
  }))
  
  # format into a character vector, one string per row; omit the df header
  specs_printed <- capture.output(print(specs_df, right=FALSE))[-1]

  # print the specs header and df rows
  cat(header)
  cat(paste0(prefix, specs_printed, "\n"), sep='')
  
  invisible(object)
}