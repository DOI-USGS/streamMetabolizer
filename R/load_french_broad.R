#' Load a short dataset from the French Broad River
#' 
#' @importFrom unitted u
#' @param attach.units logical. Should units be attached to the data.frame?
#' @return a data.frame, unitted if attach.units==TRUE
load_french_broad <- function(attach.units=TRUE) {
  file.name <- system.file("extdata", "french_broad_data.csv", package="streamMetabolizer")
  french_broad_data <- read.csv(file.name, stringsAsFactors=FALSE, skip=2)
  french_broad_units <- sapply(read.csv(file.name, stringsAsFactors=FALSE, nrow=1), function(col) as.character(col))
  french_broad_data$date.time <- as.POSIXct(strptime(french_broad_data$date.time, format="%Y-%m-%d %H:%M:%S"))
  if(attach.units) u(french_broad_data, french_broad_units) else french_broad_data
}