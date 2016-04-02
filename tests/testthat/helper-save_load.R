# saveRDS test script from http://rpubs.com/hadley/saveRDS
library(dplyr)
library(microbenchmark)
roundtrip <- function(dat, con_fun, reps=10, ...) {
  test <- tempfile()
  con <- con_fun(test, ...)
  on.exit(close(con))
  
  save <- summary(microbenchmark(saveRDS(dat, con), times=reps), unit='ms')$mean
  load <- summary(microbenchmark(x <- readRDS(test)), unit='ms')$mean
  size <- file.info(test)$size / (1024) ^ 2
  
  file.remove(test)
  data_frame(save, load, size)
}
save_load_timing <- function(dat, reps=10, ...) {
  bind_rows(
    data.frame(type='raw', level=0, roundtrip(dat, file), stringsAsFactors=FALSE),
    data.frame(type='gz', level=1, roundtrip(dat, gzfile, reps=reps, compression = 1), stringsAsFactors=FALSE),
    data.frame(type='gz', level=6, roundtrip(dat, gzfile, reps=reps, compression = 6), stringsAsFactors=FALSE),
    data.frame(type='gz', level=9, roundtrip(dat, gzfile, reps=reps, compression = 9), stringsAsFactors=FALSE),
    data.frame(type='bz', level=1, roundtrip(dat, bzfile, reps=reps, compression = 1), stringsAsFactors=FALSE),
    data.frame(type='bz', level=6, roundtrip(dat, bzfile, reps=reps, compression = 6), stringsAsFactors=FALSE),
    data.frame(type='bz', level=9, roundtrip(dat, bzfile, reps=reps, compression = 9), stringsAsFactors=FALSE),
    data.frame(type='xz', level=1, roundtrip(dat, xzfile, reps=reps, compression = 1), stringsAsFactors=FALSE),
    data.frame(type='xz', level=6, roundtrip(dat, xzfile, reps=reps, compression = 6), stringsAsFactors=FALSE),
    data.frame(type='xz', level=9, roundtrip(dat, xzfile, reps=reps, compression = 9), stringsAsFactors=FALSE)
  ) %>% 
    mutate(
      total=save+load, 
      typelevel=paste0(type, level),
      timesize=(total/max(total)) + (size/max(size))) %>%
    arrange(timesize) %>%
    mutate(typelevel=ordered(typelevel,typelevel))
}
plot_save_load_timing <- function(times) {
  ggplot(times, aes(x=typelevel, group=1)) + 
    geom_line(aes(y=timesize), color='purple', size=2) + 
    geom_line(aes(y=total*max(timesize)/max(total)), color='red') + 
    geom_line(aes(y=size*max(timesize)/max(size)), color='blue')
}
test_save_load_recovery <- function(mm) {
  # specific to metab models and uses gz6
  test <- tempfile()
  con <- gzfile(test, compression=6)
  on.exit(close(con))
  saveRDS(mm, file=con)
  mm2 <- readRDS(test)
  expect_equal(mm, mm2)
  expect_equal(get_fit(mm), get_fit(mm2))
  expect_equal(get_data(mm), get_data(mm2))
  list(original=mm, reloaded=mm2)
}