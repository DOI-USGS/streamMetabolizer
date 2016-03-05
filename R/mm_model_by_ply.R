#' Split and label data into >=24-hr days for fitting daily metabolism
#' 
#' Splits up to two data.frames, data and data_daily, into date-specific chunks.
#' These are passed to model_fun, and the results of model_fun (which must be a 
#' data.frame) are modified to include the data as a first column, then 
#' row-bound together into a single data.frame containing results from all days.
#' 
#' @param model_fun the function to apply to each data ply. This function should
#'   accept the arguments \code{c(data, data_daily, ..., day_start, day_end, 
#'   ply_date)} where \code{data_daily} is \code{NULL} when the 
#'   \code{data_daily} argument to \code{mm_model_by_ply} is missing or 
#'   \code{NULL}
#' @param data required. the data.frame containing all relevant, validated 
#'   modeling data. The solar.time column must be present.
#' @param data_daily optional. a data.frame containing inputs with a daily 
#'   timestep. The date column must be present.
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
#' @inheritParams mm_model_by_ply_prototype
#' @param ... other args to be passed through mm_model_by_ply to model_fun
#' @return a data.frame of model results
#' @import dplyr
mm_model_by_ply <- function(model_fun, data, data_daily, day_start, day_end, ...) {
  
  # Identify the data plys that will let us use a user-specified-hr window for 
  # each date (day_start to day_end, which may be != 24). store this labeling in
  # two additional columns (odd.- and even.- date.groups)
  data.plys <- as.data.frame(v(data))
  if(any(is.na(data.plys$solar.time))) stop("no values in solar.time may be NA")
  data.plys$date <- format(data.plys$solar.time, "%Y-%m-%d")
  data.plys$hour <- 24*(convert_date_to_doyhr(data.plys$solar.time) %% 1)
  
  unique.dates <- unique(as.character(data.plys$date))
  odd.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 1)]
  even.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 0)]
  
  # subtract a tiny fudge factor (one second) to the start and end bounds so we
  # don't exclude times that are essentially equal to day_start, and so we do
  # exclude times that are essentially equal to day_end
  day_start <- day_start - 1/(60*60*24)
  day_end <- day_end - 1/(60*60*24)
  
  # Assign each data row to one or two date plys. Because date plys can overlap,
  # there are two columns so that one row can belong to two dates if needed.
  dt.NA <- as.character(NA)
  for(dt in seq_along(unique.dates)) {
    dt.today <- unique.dates[dt]
    dt.yesterday <- if(dt>1) unique.dates[dt-1] else dt.NA
    dt.tomorrow <- if(dt<length(unique.dates)) unique.dates[dt+1] else dt.NA
    hr <- data.plys[data.plys$date == dt.today, 'hour']
    primary.date <- c(dt.today, dt.NA)[ifelse(hr >= day_start & hr < day_end, 1, 2)]
    secondary.date <- c(dt.yesterday, dt.NA, dt.tomorrow)[ifelse(hr <= (day_end - 24), 1, ifelse(hr < (24 + day_start), 2, 3))]
    data.plys[data.plys$date == dt.today, 'odd.date.group'] <- if(dt.today %in% odd.dates) primary.date else secondary.date
    data.plys[data.plys$date == dt.today, 'even.date.group'] <- if(dt.today %in% even.dates) primary.date else secondary.date
  } 
  
  # Apply model_fun to each ply of the data, using if-else in the lapply loop to
  # cover both the odd and even groupings in date order
  out_list <- lapply(sort(unique.dates), function(dt) {
    ply <- if(dt %in% odd.dates) {
      which(data.plys$odd.date.group == dt)
    } else {
      which(data.plys$even.date.group == dt)
    }
    out <- model_fun(
      data_ply=data.plys[ply,!(names(data.plys) %in% c('date','hour','odd.date.group','even.date.group'))], 
      data_daily_ply=if(!is.null(data_daily)) data_daily[data_daily$date == dt,] else NULL,
      day_start=day_start, day_end=day_end, ply_date=as.Date(dt, tz=tz(dt)),
      ...)
    if(is.null(out) || nrow(out)==0) NULL else data.frame(date=as.Date(dt, tz=tz(dt)), out)
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
