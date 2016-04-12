#' Split and label data into >=24-hr days for fitting daily metabolism
#' 
#' Splits up to two data.frames, data and data_daily, into date-specific chunks.
#' These are passed to model_fun. If \code{day_tests} is not empty, those 
#' validity checks are run and the results are also passed to model_fun (in 
#' \code{validity}). The results of model_fun (which must be a data.frame) are 
#' modified to include the data as a first column, then row-bound together into 
#' a single data.frame containing results from all days.
#' 
#' @param model_fun the function to apply to each data ply. This function should
#'   accept the arguments \code{c(data, data_daily, ..., day_start, day_end, 
#'   ply_date)} where \code{data_daily} is \code{NULL} when the 
#'   \code{data_daily} argument to \code{mm_model_by_ply} is missing or 
#'   \code{NULL}
#' @param data required. A data.frame to split into chunks by date, where a 
#'   'date' begins on the hour day_start and ends at the hour day_end. The 
#'   solar.time column must be present.
#' @param data_daily optional. A data.frame containing inputs with a daily 
#'   timestep, each row of which will be passed to the corresponding date chunk 
#'   from \code{data}. The date column must be present.
#' @param day_start start time (inclusive) of a day's data in number of hours 
#'   from the midnight that begins the date. For example, day_start=-1.5 
#'   indicates that data describing 2006-06-26 begin at 2006-06-25 22:30, or at 
#'   the first observation time that occurs after that time if day_start doesn't
#'   fall exactly on an observation time. For metabolism models working with 
#'   single days of input data, it is conventional/useful to begin the day the 
#'   evening before, e.g., -1.5, and to end just before the next sunrise, e.g., 
#'   30. For multiple consecutive days, it may make the most sense to start just
#'   before sunrise (e.g., 4) and to end 24 hours later. For nighttime 
#'   regression, the date assigned to a chunk of data should be the date whose 
#'   evening contains the data. The default is therefore 12 to 36 for 
#'   metab_night, of which the times of darkness will be used.
#' @param day_end end time (exclusive) of a day's data in number of hours from 
#'   the midnight that begins the date. For example, day_end=30 indicates that 
#'   data describing 2006-06-26 end at the last observation time that occurs 
#'   before 2006-06-27 06:00. See day_start for recommended start and end times.
#' @param day_tests list of tests to conduct to determine whether each date 
#'   worth of data is valid for modeling. the results of these tests will be 
#'   passed to \code{model_fun} as the \code{ply_validity} argument to that 
#'   function.
#' @param timestep_days TRUE if you would like the mean timestep length to be 
#'   calculated for each data ply and passed to \code{model_fun} as the
#'   \code{timestep_days} argument to that function. Alternatively, this may be
#'   numeric as a specifically expected timestep length in days; for example, a
#'   1-hour timestep is 1/24 is 0.0416667.
#' @param ... other args to be passed through mm_model_by_ply to model_fun
#' @return a data.frame of model results
#' @import dplyr
#' @examples
#' dat <- data_metab('10')
#' mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=2, day_end=28)$date
#' mm_model_by_ply(function(...) { data.frame(text='hi')}, 
#'   data=dat, day_start=2, day_end=28)
#' @export
mm_model_by_ply <- function(
  model_fun, data, data_daily=NULL, day_start, day_end, 
  day_tests=c('full_day', 'even_timesteps', 'complete_data'), timestep_days=TRUE, ...
) {
  
  # avoid some ugly edge cases
  if(missing(day_start) || is.null(day_start)) stop('day_start must be specified')
  if(missing(day_end) || is.null(day_end)) stop('day_end must be specified')
  if(day_end - day_start > 48) stop("day_end - day_start must not be > 48") # would break our odd/even algorithm
  if(-24 >= day_start || day_start >= 24) stop("day_start must be in (-24,24)")
  if(0 >= day_end || day_end >= 48) stop("day_end must be in (0,48)")
  
  # Identify the data plys that will let us use a user-specified-hr window for 
  # each date (day_start to day_end, which may be != 24). store this labeling in
  # two additional columns (odd.- and even.- date.groups)
  data.plys <- as.data.frame(v(data))
  if(!('solar.time' %in% names(data.plys))) stop("data must contain a 'solar.time' column")
  if(any(is.na(data.plys$solar.time))) stop("no values in solar.time may be NA")
  if((min_timestep <- mm_get_timestep(data$solar.time, format='unique')[1]) <= 0) 
    stop("min timestep is <= 0: ", min_timestep, " days")
  if(!is.null(data_daily)) {
    if(!('date' %in% names(data_daily))) stop("data_daily must contain a 'date' column")
    min_datestep <- mm_get_timestep(data_daily$date, format='unique')
    if(length(min_datestep) > 0 && min_datestep[1] <= 0)
      stop("min datestep is <= 0: ", min_datestep, " days")
  }
  doyhr <- convert_date_to_doyhr(data.plys$solar.time)
  data.plys$date <- as.character(as.Date(lubridate::floor_date(data.plys$solar.time, 'year')) + as.difftime(floor(doyhr)-1, units='days'))
  data.plys$hour <- 24*(doyhr %% 1)
  
  # subtract a tiny fudge factor (one second) to the start and end bounds so we
  # don't exclude times that are essentially equal to day_start, and so we do
  # exclude times that are essentially equal to day_end
  day_start_fudge <- day_start - 1/(60*60*24)
  day_end_fudge <- day_end - 1/(60*60*24)
  
  # This function speeds things up in both of the next two loops. 'runs' refers
  # to sequential runs of values in date.group that have the same value.
  get_runs <- function(date.group) {
    values <- '.dplyr.var'
    replace(date.group, is.na(date.group), '') %>%
      rle() %>%
      unclass %>% 
      as_data_frame %>% 
      mutate(end=cumsum(lengths), start=c(1, 1+end[-n()])) %>%
      filter(values != '') %>%
      select(date=values, start, end)
  }
  # Determine which rows go with which date, and label dates as either odd or
  # even (no relation to the actual value of the date)
  date_rows <- get_runs(data.plys$date)
  odd.dates <- date_rows[which(seq_len(nrow(date_rows)) %% 2 == 1),]$date
  even.dates <- date_rows[which(seq_len(nrow(date_rows)) %% 2 == 0),]$date
  
  # Assign each data row to one or two date plys. Because date plys can overlap,
  # there are two columns so that one row can belong to two dates if needed.
  dt.NA <- as.character(NA)
  data.plys <- bind_cols(
    data.plys,
    bind_rows(lapply(seq_len(nrow(date_rows)), function(dt) {
      run <- date_rows[dt,]
      dt.today <- run$date
      dt.yesterday <- if(dt>1) date_rows[[dt-1,'date']] else dt.NA
      dt.tomorrow <- if(dt<nrow(date_rows)) date_rows[[dt+1,'date']] else dt.NA
      data.plys.rows <- run$start:run$end
      hr <- data.plys[data.plys.rows, 'hour']
      primary.date <- c(dt.today, dt.NA)[ifelse(hr >= day_start_fudge & hr < day_end_fudge, 1, 2)]
      secondary.date <- c(dt.yesterday, dt.NA, dt.tomorrow)[ifelse(hr <= (day_end_fudge - 24), 1, ifelse(hr < (24 + day_start_fudge), 2, 3))]
      if(dt.today %in% odd.dates) {
        data_frame(odd.date.group=primary.date, even.date.group=secondary.date)
        #data.plys[data.plys.rows, c('odd.date.group','even.date.group')] <- c(primary.date, secondary.date)
      } else { # dt.today %in% even.dates
        data_frame(odd.date.group=secondary.date, even.date.group=primary.date)
        #data.plys[data.plys.rows, c('odd.date.group','even.date.group')] <- c(secondary.date, primary.date)
      }
    }))) %>%
    as.data.frame(stringsAsFactors=FALSE)
  
  # filter out dates that don't ever appear on their own - this especially weeds
  # out first and last dates that were probably meant just to be part of the
  # second or penultimate dates when the date range is >24
  date.pairings <- setNames(unique(data.plys[, c('odd.date.group','even.date.group')]), c('a','b'))
  date.pairings <- rbind(date.pairings, setNames(date.pairings, c('b','a')))
  date.pairings <- date.pairings[!is.na(date.pairings$a),]
  solo.dates <- date.pairings[is.na(date.pairings$b),]$a
  tbl.dates <- table(date.pairings$a)
  double.dates <- names(tbl.dates)[tbl.dates > 1]
  unique.dates <- sort(unique(c(solo.dates, double.dates)))
  
  # Apply model_fun to each ply of the data, using if-else in the lapply loop to
  # cover both the odd and even groupings in date order
  runs <- bind_rows(
    get_runs(data.plys$odd.date.group),
    get_runs(data.plys$even.date.group)) %>%
    filter(date %in% unique.dates) %>%
    arrange(date) %>%
    mutate(date=as.Date(date, tz=tz(data$solar.time)))
  hour <- odd.date.group <- even.date.group <- '.dplyr.var'
  data_plys <- data.plys %>%
    select(-date, -hour, -odd.date.group, -even.date.group)
  out_list <- lapply(seq_along(runs$date), function(dt) {
    # pick out the inst & daily plys for this date
    run <- runs[dt,]
    data_ply <- data_plys[run$start:run$end,]
    data_daily_ply <- 
      if(!is.null(data_daily)) {
        data_daily[match(run$date, data_daily$date),]
      } else {
        NULL
      }
    ply_date <- run$date
    
    # compute timestep and run validity checks if requested. if timestep_days is
    # FALSE or NA, it will be passed as NA to the model_fun and only computed
    # there if needed for specific tests
    if(length(timestep_days) > 1) stop("expecting no more than 1 value in timestep_days")
    timestep_days <- if(isTRUE(timestep_days)) {
      mm_get_timestep(data_ply$solar.time, format='mean') 
    } else if(is.na(timestep_days) || timestep_days==FALSE) {
      NA
    } else timestep_days
    ply_validity <- if(length(day_tests) > 0) {
      mm_is_valid_day(
        data_ply=data_ply, day_start=day_start, day_end=day_end, day_tests=day_tests, 
        ply_date=ply_date, timestep_days=timestep_days)
    } else NA
    
    # run the user's model_fun
    out <- model_fun(
      data_ply=data_ply, data_daily_ply=data_daily_ply,
      day_start=day_start, day_end=day_end, ply_date=ply_date,
      ply_validity=ply_validity, timestep_days=timestep_days, 
      ...)
    # attach a date column if anything was returned
    if(is.null(out) || nrow(out)==0) NULL else data.frame(date=ply_date, out)
  })
  # combine the 1-row dfs, making sure to start with a 0-row df that represents 
  # the dfs with the most columns to keep the column ordering consistent across 
  # runs. dplyr::bind_rows should take care of any missing columns, since these
  # are matched by name, and missing columns are filled with NAs
  count_cols <- function(out) { if(is.null(out)) 0 else ncol(out) }
  example_choice <- which.max(sapply(out_list, count_cols))
  if(length(example_choice) == 0) {
    data.frame(date=as.Date(NA)) 
  } else {
    out_example <- out_list[[example_choice]][FALSE,]
    bind_rows(c(list(out_example), out_list)) %>% as.data.frame() 
  }
}
