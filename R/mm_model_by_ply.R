#' Split and label data into >=24-hr days for fitting daily metabolism
#' 
#' @param data the data.frame containing all relevant, validated modeling data
#' @param model_fun the function to apply to each data ply
#' @inheritParams mm_is_valid_day
#' @param ... additional args passed to model_fun (data, day_start, and day_end
#'   will also be passed)
#' @import dplyr
#' @return a data.frame of fitting results
mm_model_by_ply <- function(data, model_fun, day_start, day_end, ...) {
  
  # my use of dplyr might be responsible for seg faults here. try it without
  # dplyr.
  use.dplyr <- FALSE
  
  # Identify the data plys that will let us use a user-specified-hr window for
  # each date (day_start to day_end, which may be != 24) - this labeling can
  # be stored in two additional columns (odd.- and even.- date.groups)
  if(use.dplyr) {
    local.time <- hour <- ".dplyr.var"
    data.plys <- data %>% 
      mutate(date=as.Date(format(local.time, "%Y-%m-%d")),
             hour=24*(convert_date_to_doyhr(local.time) %% 1))
  } else {
    data.plys <- as.data.frame(v(data))
    data.plys$date <- as.Date(format(data.plys$local.time, "%Y-%m-%d"))
    data.plys$hour <- 24*(convert_date_to_doyhr(data.plys$local.time) %% 1)
  }
  
  unique.dates <- unique(data.plys$date)
  odd.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 1)]
  even.dates <- unique.dates[which(seq_along(unique.dates) %% 2 == 0)]
  
  if(use.dplyr) {
    data.plys <- data.plys %>% 
      group_by(date) %>%
      mutate(odd.date.group=if(date[1] %in% odd.dates) date else c(date[1]-1, as.Date(NA), date[1]+1)[ifelse(hour <= (day_end - 24), 1, ifelse(hour < (24 + day_start), 2, 3))],
             even.date.group=if(date[1] %in% even.dates) date else c(date[1]-1, as.Date(NA), date[1]+1)[ifelse(hour <= (day_end - 24), 1, ifelse(hour < (24 + day_start), 2, 3))]) %>%
      ungroup() %>% select(-date)
  } else {
    for(dt in unique(data.plys$date)) {
      hr <- data.plys[data.plys$date == dt, 'hour']
      primary.date <- c(dt, as.Date(NA))[ifelse(hr > day_start & hr <= day_end, 1, 2)]
      secondary.date <- c(dt-1, as.Date(NA), dt+1)[ifelse(hr <= (day_end - 24), 1, ifelse(hr < (24 + day_start), 2, 3))]
      data.plys[data.plys$date == dt, 'odd.date.group'] <- if(dt %in% odd.dates) primary.date else secondary.date
      data.plys[data.plys$date == dt, 'even.date.group'] <- if(dt %in% even.dates) primary.date else secondary.date
    } 
    data.plys <- data.plys[,-which(names(data.plys)=='date')]
  }
  
  # Estimate daily metabolism for each ply of the data, using two group_by/do
  # combinations to cover the odd and even groupings
  if(use.dplyr) {
    . <- odd.date.group <- even.date.group <- ".dplyr.var"
    out.all <- 
      bind_rows(
        data.plys %>% group_by(date=odd.date.group) %>% do(model_fun(., ..., day_start=day_start, day_end=day_end)), # filter(!is.na(odd.date.group)) %>%, filter(!is.na(even.date.group)) %>% 
        data.plys %>% group_by(date=even.date.group) %>% do(model_fun(., ..., day_start=day_start, day_end=day_end))) %>% 
      filter(!is.na(date), date %in% unique.dates) %>%
      arrange(date) 
  } else {
    out.all <- do.call(
      rbind, 
      lapply(sort(unique.dates), function(dt) {
        ply <- if(dt %in% odd.dates) {
          which(data.plys$odd.date.group == dt)
        } else {
          which(data.plys$even.date.group == dt)
        }
        out <- model_fun(data.plys[ply,!(names(data.plys) %in% c('odd.date.group','even.date.group'))], 
                         ..., day_start=day_start, day_end=day_end)
        out$date <- dt
        out[,c(ncol(out), 1:(ncol(out)-1))]
      }))
  } 
  
  out.all
}
