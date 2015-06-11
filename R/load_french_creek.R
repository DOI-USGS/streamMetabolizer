#' Load a short dataset from French Creek
#' 
#' @importFrom unitted u
#' @param attach.units logical. Should units be attached to the data.frame?
#' @return a data.frame, unitted if attach.units==TRUE
load_french_creek <- function(attach.units=TRUE) {
  file.name <- system.file("extdata", "french_creek_data.csv", package="streamMetabolizer")
  french_creek_data <- read.csv(file.name, stringsAsFactors=FALSE, skip=2)
  french_creek_units <- sapply(read.csv(file.name, stringsAsFactors=FALSE, nrow=1), function(col) as.character(col))
  french_creek_data$date.time <- as.POSIXct(strptime(french_creek_data$date.time, format="%m/%d/%Y %H:%M"), tz="GMT")
  if(attach.units) u(french_creek_data, french_creek_units) else french_creek_data
}