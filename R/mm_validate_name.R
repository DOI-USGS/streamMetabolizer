#' Check the validity of a model name
#' 
#' Check the syntactic & scientific validity of a model name. Returns the model
#' name if it's valid, otherwise gives an error
#' 
#' @inheritParams specs
#' @examples 
#' mm_validate_name("b_np_oipi_tr_plrckm.stan")
#' \dontrun{
#' mm_validate_name("b_np_oipn") # throws error
#' }
#' @export
mm_validate_name <- function(model_name) {
  
  # under lazy eval, model_name might not break until here (e.g., 
  # mm_validate_name(mm_name('mle', badarg=4))). force it here to give mm_name a
  # chance to throw its own errors.
  model_name <- force(model_name)
  
  # require parseable name
  parse_problem <- c()
  parsed <- tryCatch(
    mm_parse_name(model_name),
    error=function(e) { 
      parse_problem <<- "could not parse model name. try constructing with mm_name()"
    })
  if(length(parse_problem) > 0) {
    stop(parse_problem)
  }
  
  # require valid type
  type <- parsed$type
  valid_types <- eval(formals(mm_name)$type)
  if(is.na(type) || !(type %in% valid_types)) {
    stop('model name implies unknown model type (', type, '). try constructing with mm_name()')
  }
  
  # check against known or findable valid names
  valid_names <- mm_valid_names(type)
  if(basename(model_name) != model_name) {
    mm_locate_filename(model_name)
  } else if(!(model_name %in% valid_names)) {
    stop("model_name (", model_name, ") is not among valid ", type, 
         sprintf(" model_names (see mm_valid_names('%s'))", type))
  }
  
  # return the model name
  model_name
}
