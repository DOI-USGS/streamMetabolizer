context("output of metabolism models")

# use a subset of data from Bob
library(dplyr)
library(unitted)
french <- streamMetabolizer:::load_french_creek()

test_that("metabolism models run & produce reasonable output", {
  
  # metab_mle
  mm <- metab_mle(data=v(french))
  expect_is(mm, "metab_mle")
  expect_is(slot(mm, "fit"), "data.frame") # specific to this model
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")

  mm <- metab_night(data=v(french))
  
  mm <- metab_mle(data=v(french))
  
  dontrun <- function() {
    library(dplyr)
    library(unitted)
    library(mda.streams)
    config <- read.table("../stream_metab_usa/p2_metab/out/150714 0.0.2 local_makefile_run/condor_config.tsv", header=T, sep="\t", colClasses="character")
    data <- config_to_data(config[40,], 40, metab_night, list())
    data_ply <- data[which(as.Date(v(data$local.time)+as.difftime(12,units="hours")) %in% as.Date(paste0("2014-04-",c("10","11","12")))),]
    data_ply <- u(data_ply, get_units(data_ply) %>% replace(., which(.=="mg L^-1"), "mgO2 L^-1"))
    mm <- metab_night(v(data_ply))
    mm <- metab_night(data_ply)
  }
  
})

test_that("metabolism predictions (predict_metab, predict_DO) make sense", {
  
  # metab_mle
  mm <- metab_bayes(data=v(french), adaptSteps=100, burnInSteps=400, numSavedSteps=4000)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # library(ggplot2); ggplot(DO_preds, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  # library(ggplot2); ggplot(DO_preds_Aug24, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  
})

test_that("metabolism predictions (predict_metab, predict_DO) make sense", {
  
  # metab_mle
  mm <- metab_mle(data=v(french))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # library(ggplot2); ggplot(DO_preds, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  # library(ggplot2); ggplot(DO_preds_Aug24, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  
})

test_that("metabolism models can be fit with K specified", {
  
  # metab_mle
  K600 <- data.frame(date=unique(as.Date(v(french)$local.time)), K600=c(NA, 30, 30, 0, NA, 0, 50, 40))
  mm <- metab_mle(data=v(french), K600=K600)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # library(ggplot2); ggplot(DO_preds, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  # library(ggplot2); ggplot(DO_preds_Aug24, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  
})


test_that("metab_models can be saved", {
  
  # save and reload metab_mle
  mm <- metab_mle(data=v(french))
  save(mm, file="test_temp.RData")
  mm_orig <- mm
  load("test_temp.RData")
  file.remove("test_temp.RData")
  expect_equal(mm_orig, mm)
  expect_equal(get_fit(mm_orig), get_fit(mm))
  
})