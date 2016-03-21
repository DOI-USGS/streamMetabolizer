#' Get a demo dataset for modeling metabolism
#' 
#' Get a formatted data.frame of inputs from which metabolism can be modeled. 
#' These test data were provided by Bob Hall.
#' 
#' @param num_days the number of days to include in the data. character format 
#'   because only certain numbers of days are permitted (see defaults in Usage 
#'   for the accepted options)
#' @param res character specifying the desired resolution of the data in minutes
#'   (character; see defaults in Usage for the accepted options)
#' @param flaws character specifying one or more flaws to include in the data, or
#'   empty (\code{c()}) for no flaws. default is no flaws.
#' @inheritParams mm_model_by_ply
#' @inheritParams load_french_creek
#' @importFrom unitted u v get_units
#' @examples 
#' head(data_metab())
#' head(data_metab(res='30'))
#' @export
data_metab <- function(
  num_days=c('1','3','10'), 
  res=c('5','10','15','30'),
  flaws=c('missing middle', 'missing start', 'missing end', 'missorted', 'duplicated'),
  day_start=4, day_end=28, attach.units=FALSE) {
  
  # check inputs
  num_days <- match.arg(num_days)
  res <- as.numeric(match.arg(res))
  flaws <- if(missing(flaws)) c() else match.arg(flaws, several.ok=TRUE)
  
  # start with the same data every time
  french <- load_french_creek(attach.units=attach.units)
  # take off units temporarily to make it easier to manipulate datetimes
  if(attach.units) {
    french_units <- get_units(french)
    french <- v(french)
  }
  french <- french[order(french$solar.time),]
  # fill in holes in the part of the data we'll be using
  french <- french[c(1:6352, rep(6353, 3), 6354:7772, rep(7773, 2), 7774:nrow(french)),] 
  french[6354:6355,'solar.time'] <- french[6353,'solar.time'] + as.difftime(c(5,10), units='mins')
  french[7776,'solar.time'] <- french[7775,'solar.time'] + as.difftime(5, units='mins')
  # save dates vec now to help with res changes later
  orig_times <- french$solar.time 
  
  # subset by num_days
  date_start <- "2012-09-18"
  date_end <- format(as.Date(date_start) + as.numeric(num_days)-1, "%Y-%m-%d")
  french <- mm_filter_dates(french, date_start=date_start, date_end=date_end, day_start=day_start, day_end=day_end)
  
  # add flaws
  day_length <- 24 * 12 # 12 obs/hr (every 5 mins) in raw data
  day2_start <- 1 + day_length
  day3_start <- day2_start + day_length
  skip_rows <- c()
  if('missing start' %in% flaws) {
    skip_rows <- c(
      skip_rows, 
      switch(
        num_days, 
        '1' = 1:6, 
        '3' = day2_start + 0:23, 
        '10' = c(day3_start + 0:41, (day3_start + 2*day_length) + 0:17)))
  }
  if('missing middle' %in% flaws) {
    skip_rows <- c(
      skip_rows, 
      switch(
        num_days, 
        '1' = 31:48, 
        '3' = day2_start + 31:60, 
        '10' = c(day3_start + 71:130, (day3_start + 2*day_length) + 14:55)))
  }
  if('missing end' %in% flaws) {
    skip_rows <- c(
      skip_rows, 
      switch(
        num_days, 
        '1' = day2_start - 1:6, 
        '3' = day3_start + 1:24, 
        '10' = day_length + c(day3_start - 1:42, (day3_start + 2*day_length) - 1:18)))
  }
  if(length(skip_rows) > 0) 
    french <- french[-skip_rows, ]

  if('missorted' %in% flaws) {
    swap_rows <- switch(
      num_days,
      '1' = list(from=31:48, to=66),
      '3' = list(from=day2_start + 31:60, to=day2_start + 42),
      '10' = list(from=c(day3_start + 71:130, (day3_start + 2*day_length) + 14:55), to=day2_start + 45))
    start_rows <- 1:(swap_rows$to - 1)
    middle_rows <- swap_rows$from
    end_rows <- swap_rows$to:nrow(french)
    french <- french[c(start_rows, middle_rows, end_rows), ]
  }
  
  if('duplicated' %in% flaws) {
    dup_rows <- switch(
      num_days, 
      '1' = 31:48, 
      '3' = day2_start + 31:60, 
      '10' = c(day3_start + 71:130, (day3_start + 2*day_length) + 14:55))
    french <- french[sort(c(1:nrow(french), dup_rows)), ]
  }

  # subset by resolution (yes, before flaws=duplicated)
  sub_times <- orig_times[seq(1, length(orig_times), by=res/5)]
  french <- french[french$solar.time %in% sub_times, ]

  # add back units if requested
  if(attach.units) {
    french <- u(french, french_units)
  }
  
  # return
  french
}
