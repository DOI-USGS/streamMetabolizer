#' Return the average timestep in days
#' 
#' @param datetimes a vector of date-times in POSIXct format from which to 
#'   compute the average timestep
#' @param format the format in which to return the timestep. 'mean' always 
#'   returns one value; 'unique' may return more than one depending on the 
#'   variation in timesteps and the value of \code{digits}.
#' @param require_unique logical. should it be required that there is exactly
#'   one unique timestep (within the given tolerance \code{tol})?
#' @param tol if \code{format == 'unique'}, unique values are first calculated 
#'   to machine precision, but then subsetted to those that differ from one 
#'   another by at least tol, where tol is a time difference in units of days 
#'   (and thus 1/(24*60*60) is one second).
#' @importFrom unitted v
#' @examples {
#' datetimes <- Sys.time()+ as.difftime(c(0,304,600,900.2,1200,1500,1800), units='secs')
#' mm_get_timestep(datetimes, 'unique', tol=1/(24*60*60))
#' mm_get_timestep(datetimes, 'unique', tol=5/(24*60*60))
#' mm_get_timestep(datetimes, 'unique', tol=10/(24*60*60))
#' mm_get_timestep(datetimes, 'unique', tol=300/(24*60*60))
#' mm_get_timestep(datetimes, 'mean')
#' mm_get_timestep(datetimes, 'mean', require_unique=TRUE, tol=300/(24*60*60))
#' datetimes <- Sys.time()+ as.difftime(c(-1,0,2,4,5,6,8,10), units='days')
#' mm_get_timestep(datetimes, 'modal')
#' mm_get_timestep(c(), 'mean')
#' mm_get_timestep(c(), 'unique')
#' mm_get_timestep(c(), 'modal')
#' \dontrun{
#' # all of these should and do give errors:
#' mm_get_timestep(datetimes, 'mean', require_unique=TRUE, tol=1/(24*60*60))
#' mm_get_timestep(datetimes, 'unique', tol=5/(24*60*60), require_unique=TRUE)
#' mm_get_timestep(c(), 'mean', require_unique=TRUE)
#' mm_get_timestep(c(), 'unique', require_unique=TRUE)
#' mm_get_timestep(c(), 'modal', require_unique=TRUE)
#' }
#' }
#' @export
mm_get_timestep <- function(datetimes, format=c('mean','unique','modal'), require_unique=FALSE, tol=60/(24*60*60)) {
  timesteps <- as.numeric(diff(v(datetimes)), units="days")
  timestep <- switch(
    match.arg(format),
    mean = {
      if(require_unique == TRUE)
        mm_get_timestep(datetimes, format='unique', require_unique=TRUE, tol=tol)
      if(length(timesteps) == 0) 
        c() 
      else 
        mean(timesteps, na.rm=TRUE)
    },
    unique = {
      all_unique <- sort(unique(timesteps))
      sufficiently_unique <- c()
      while(length(all_unique) > 0) {
        sufficiently_unique <- c(sufficiently_unique, all_unique[1])
        all_unique <- all_unique[which(all_unique > tail(sufficiently_unique, 1) + tol)]
      }
      if(require_unique == TRUE && length(sufficiently_unique) != 1)
        stop('!=1 unique timestep')
      sufficiently_unique
    },
    modal = {
      all_unique <- sort(unique(timesteps))
      tbl_unique <- sapply(all_unique, function(tstep) sum(timesteps == tstep), USE.NAMES=FALSE)
      sufficiently_unique <- mm_get_timestep(
        datetimes=cumsum(c(0,all_unique)), format='unique', require_unique=require_unique, tol=tol)
      if(length(sufficiently_unique) == 0) return(c())
      round_unique <- approx(sufficiently_unique, sufficiently_unique, xout=all_unique, method='constant', rule=2)$y
      tbl_suff_unique <- sapply(sufficiently_unique, function(ru) sum(tbl_unique[which(round_unique == ru)]))
      all_suff_unique <- sapply(sufficiently_unique, function(ru) {
        matches <- which(round_unique == ru)
        best_match <- matches[which.max(tbl_unique[matches])]
        all_unique[best_match]
      })
      all_suff_unique[which.max(tbl_suff_unique)]
    })
  timestep
}
