#' Split and label data into >=24-hr days for fitting daily metabolism
#' 
#' Splits up to two data.frames, data and data_daily, into date-specific chunks.
#' These are passed to model_fun, and the results of model_fun (which must be a 
#' data.frame) are modified to include the data as a first column, then 
#' row-bound together into a single data.frame containing results from all days.
#' 
#' @param data required. the data.frame containing all relevant, validated 
#'   modeling data. The local.time column must be present.
#' @param data_daily optional. a data.frame containing inputs with a daily 
#'   timestep. The local.date column must be present.
#' @param model_fun the function to apply to each data ply. This function should
#'   accept the arguments \code{c(data, data_daily, ..., day_start, day_end, local_date)}
#'   where \code{data_daily} is \code{NULL} when the \code{data_daily} argument
#'   to \code{mm_model_by_ply} is missing or \code{NULL}
#' @inheritParams mm_is_valid_day
#' @param ... additional args passed to model_fun (data, day_start, and day_end 
#'   will also be passed)
#' @return a data.frame of fitting results
mm_model_by_ply <- function(data, data_daily, model_fun, day_start, day_end, ...) {
  
  # Identify the data plys that will let us use a user-specified-hr window for 
  # each date (day_start to day_end, which may be != 24). store this labeling in
  # two additional columns (odd.- and even.- date.groups)
  data.plys <- as.data.frame(v(data))
  data.plys$local.date <- as.Date(format(data.plys$local.time, "%Y-%m-%d"))
  data.plys$hour <- 24*(convert_date_to_doyhr(data.plys$local.time) %% 1)
  
  unique.dates <- unique(data.plys$local.date)
  odd.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 1)]
  even.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 0)]
  
  # add a tiny fudge factor to the start and end bounds so we don't exclude
  # times that are essentially equal to day_start or day_end
  day_start <- day_start - 0.0001
  day_end <- day_end + 0.0001
  
  for(dt in unique(data.plys$local.date)) {
    hr <- data.plys[data.plys$local.date == dt, 'hour']
    primary.date <- c(dt, as.Date(NA))[ifelse(hr >= day_start & hr <= day_end, 1, 2)]
    secondary.date <- c(dt-1, as.Date(NA), dt+1)[ifelse(hr <= (day_end - 24), 1, ifelse(hr < (24 + day_start), 2, 3))]
    data.plys[data.plys$local.date == dt, 'odd.date.group'] <- if(dt %in% odd.dates) primary.date else secondary.date
    data.plys[data.plys$local.date == dt, 'even.date.group'] <- if(dt %in% even.dates) primary.date else secondary.date
  } 
  
  # Estimate daily metabolism for each ply of the data, using two
  # group_by/do==lapply loops to cover the odd and even groupings
  out.all <- do.call(
    rbind, 
    lapply(sort(unique.dates), function(dt) {
      ply <- if(dt %in% odd.dates) {
        which(data.plys$odd.date.group == dt)
      } else {
        which(data.plys$even.date.group == dt)
      }
      out <- model_fun(
        data=data.plys[ply,!(names(data.plys) %in% c('local.date','hour','odd.date.group','even.date.group'))], 
        data_daily=if(!missing(data_daily) && !is.null(data_daily)) data_daily[data_daily$local.date == dt,] else NULL,
        ..., day_start=day_start, day_end=day_end, local_date=dt)
      if(is.null(out) || nrow(out)==0) NULL else data.frame(local.date=dt, out)
    }))
  
  out.all
}
