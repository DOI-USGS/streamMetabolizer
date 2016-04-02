context("model_by_ply")

test_that("mm_model_by_ply creates intuitive ply_dates from day_start and day_end", {
  
  dat <- data_metab('10', res='30')
  known_dates <- mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=4, day_end=28)$date
  expect_equal(length(known_dates), 10) # 9/18 to 9/27
  
  # still return those 10 dates for non-conventional date ranges
  dat <- data_metab('10', res='30', day_start=-12, day_end=35.5)
  mmp <- mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=-12, day_end=35.5)
  expect_equal(known_dates, mmp$date) # almost 48
  expect_equal(mmp$date, as.Date(mmp$data_ply_start) + 1)
  expect_equal(mmp$date, as.Date(mmp$data_ply_end) - 1)
  mmp <- mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=4, day_end=7)
  expect_equal(known_dates, mmp$date[-11]) # bonus day b/c dat is spread wide, day_start is near day_end
  expect_equal(mmp$date, as.Date(mmp$data_ply_start))
  expect_equal(mmp$date, as.Date(mmp$data_ply_end))
  
  # still return those 10 dates for windows centered on next date
  dat <- data_metab('10', res='30', day_start=20, day_end=39)
  mmp <- mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=22, day_end=35)
  expect_equal(known_dates, mmp$date) # centered on next date
  expect_equal(mmp$date, as.Date(mmp$data_ply_start) + 0)
  expect_equal(mmp$date, as.Date(mmp$data_ply_end) - 1)
  expect_error(mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=25, day_end=30), "day_start must be in (-24,24)", fixed=TRUE) # can't be all on next date
  
  # still return those 10 dates for windows centered on previous date
  dat <- data_metab('10', res='30', day_start=-13, day_end=6)
  mmp <- mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=-12, day_end=5)
  expect_equal(known_dates, mmp$date) # centered on prev date
  expect_equal(mmp$date, as.Date(mmp$data_ply_start) + 1)
  expect_equal(mmp$date, as.Date(mmp$data_ply_end) + 0)
  expect_error(mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=-12, day_end=-2), "day_end must be in (0,48)", fixed=TRUE) # can't be all on prev date

  # catch real problems
  dat <- data_metab('10', res='30')
  expect_error(mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=-22, day_end=28), 'day_end - day_start must not be > 48')
  expect_error(mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=-26, day_end=4), 'day_start must be in (-24,24)', fixed=TRUE)
  expect_error(mm_model_by_ply(mm_model_by_ply_prototype, data=dat, day_start=22, day_end=49), 'day_end must be in (0,48)', fixed=TRUE)
})

