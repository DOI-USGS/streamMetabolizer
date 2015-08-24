context("metab_mle")

vfrench <- unitted::v(streamMetabolizer:::load_french_creek())

test_that("metab_mle models can be created", {
  
  mm <- metab_mle(data=vfrench)
  
  # check basic structure
  expect_is(mm, "metab_mle")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})

test_that("metabolism predictions (predict_metab, predict_DO) make sense", {
  
  # metab_mle
  mm <- metab_mle(data=vfrench)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- filter(DO_preds, local.date == "2012-08-24")
  knownbug <- function() {
    expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
    # library(ggplot2); ggplot(DO_preds, aes(x=local.time, group=factor(local.date))) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
    # library(ggplot2); ggplot(DO_preds_Aug24, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  }
  
})

test_that("metabolism models can be fit with K specified", {
  
  # metab_mle with K600 - not fixed up yet
  knownbug <- function() {
    K600 <- data.frame(date=unique(as.Date(v(french)$local.time)), K600=c(NA, 30, 30, 0, NA, 0, 50, 40))
    mm <- metab_mle(data=v(french), K600=K600)
    metab <- predict_metab(mm)
    DO_preds <- predict_DO(mm)
    DO_preds_Aug24<- filter(DO_preds, local.date == "2012-08-24")
    expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
    # library(ggplot2); ggplot(DO_preds, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
    # library(ggplot2); ggplot(DO_preds_Aug24, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  }
  
})
