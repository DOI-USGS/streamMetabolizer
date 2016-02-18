context("model_by_ply")

test_that("mm_model_by_ply", {
  
  ### Define dummy functions to demo & test the process ###
  
  #' Accept high-level args and call mm_model_by_ply to get things done
  dummy_outer_fun <- function(
    data, data_daily, day_start, day_end, 
    favorite # inheritParams dummy_inner_fun
  ) {
    streamMetabolizer:::mm_model_by_ply(
      dummy_middle_fun, data=data, data_daily=data_daily, # for mm_model_by_ply
      day_start=day_start, day_end=day_end, # for mm_model_by_ply
      favorite=favorite) # for dummy_inner_fun
  }
  
  #' Get called by mm_model_by_ply; call dummy_inner_fun to get answers
  dummy_middle_fun <- function(
    data_ply, data_daily_ply, day_start, day_end, solar_date, # supplied by mm_model_by_ply; inheritParams mm_model_by_ply_prototype
    favorite # inheritParams dummy_inner_fun
  ) {
    # inst & daily data should both have been subsetted to just one date (or none if there isn't a match)
    expect_less_than(nrow(data_ply), nrow(dat))
    expect_less_than(nrow(data_daily_ply), 2)
    if(nrow(data_ply) > 0 && nrow(data_daily_ply) > 0) {
      expect_true(data_daily_ply$solar.date %in% unique(as.Date(data_ply$solar.time)))
    }
    
    # fix up empty daily data because we can do that here (or in dummy_inner_fun)
    if(nrow(data_daily_ply) == 0) data_daily_ply <- data.frame(solar.date=solar_date, daily_row=NA)
      
    # call an inner function to rename arguments
    dummy_inner_fun(inst_dat=data_ply, daily_dat=data_daily_ply, favorite)
  }
  
  #' Combine the three data.frames into one and return
  dummy_inner_fun <- function(inst_dat, daily_dat, favorite) {
    
    # favorite should have been passed through untouched
    expect_equal(nrow(favorite), 10)

    # fix up empty inst data because we can do that here (or in dummy_middle_fun)
    if(nrow(inst_dat) == 0) inst_dat <- data.frame(solar.time=as.POSIXct(as.character(daily_dat$solar.date), tz='UTC'), inst_row=NA)
        
    # return combined df so we can inspect it
    data.frame(inst=inst_dat, daily=daily_dat, fav=favorite[daily_dat$daily_row,][1,])
  }
  
  ### Run the dummy functions & check output ###
  
  # define input data
  dat <- data.frame(solar.time=(as.POSIXct(as.character(Sys.Date()), tz='UTC') + as.difftime(1:(24*7), units="hours")), inst_row=1:(24*7))
  datd <- data.frame(solar.date=(Sys.Date() + as.difftime(c(-2,-1,1:3,5), units="days")), daily_row=1:6) #intentionally create date mismatch w/ dat
  extra_dailies <- (Sys.Date() + as.difftime(c(-2,-1), units="days"))
  missing_dailies <- (Sys.Date() + as.difftime(c(0,4,6,7), units="days"))
  fav <- data.frame(color=c("red","blue"), scale="b minor", digit=0:9)
  
  # run & checks: short day
  out <- dummy_outer_fun(data=dat, data_daily=datd, day_start=5, day_end=21, favorite=fav)
  expect_equal(nrow(out), 120)
  expect_equal(as.Date(out$inst.solar.time), out$daily.solar.date) # daily & inst data should be matched up right by date
  expect_equal(out$solar.date, out$daily.solar.date) # mm_model_by_ply date & daily date should be equal
  expect_equal(missing_dailies, unique(out[is.na(out$daily.daily_row),'solar.date'])) # inst obs w/o daily data should be kept
  expect_false(any(extra_dailies %in% unique(out$solar.date))) # daily data w/o inst obs should be removed
  out_hours <- as.numeric(format(out$inst.solar.time, "%H"))
  expect_equal(min(out_hours[!is.na(out$inst.inst_row)]), 5) # excepting NA inst days, minimum hour of day should be day_start=5
  expect_equal(max(out_hours), 21) # maximum hour of day should be day_end=21
  
  # run & checks: long day
  out <- dummy_outer_fun(data=dat, data_daily=datd, day_start=5, day_end=35, favorite=fav)
  expect_equal(nrow(out), 207)
  expect_true(all(as.Date(out$inst.solar.time) >= out$solar.date - as.difftime(1, units="days") & 
                    as.Date(out$inst.solar.time) <= out$solar.date + as.difftime(1, units="days"))) # were daily & inst data matched up right?
  out_hours <- as.numeric(format(out$inst.solar.time, "%H"))
  expect_true(all(table(out_hours)[as.character(5:11)] == 13)) # 5am to 11am (=35th hour) should be nearly twice as frequent as afternoon-evening
  expect_true(all(table(out_hours)[as.character(12:23)] == 7)) # afternoon-evening should be represented once per day for the 7 days in dat
  
})