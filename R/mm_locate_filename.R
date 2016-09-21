#' Look for a model file
#' 
#' Looks first in the models folder of the streamMetabolizer package, second 
#' along the relative or absolute file path given by model_name
#' 
#' @param model_name a model file in the 'models' folder of the 
#'   streamMetabolizer package or a relative or absolute file path of a model 
#'   file
#' @return a file path if the file exists or an error otherwise
#' @keywords internal
mm_locate_filename <- function(model_name) {
  
  package_path <- system.file(paste0("models/", model_name), package="streamMetabolizer")
  other_path <- model_name
  
  if(file.exists(package_path)) return(package_path)
  if(file.exists(other_path)) return(other_path)
  
  stop("could not locate the model file at ", package_path, " or ", other_path) 
}