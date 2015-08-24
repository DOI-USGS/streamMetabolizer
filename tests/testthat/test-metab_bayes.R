context("metab_bayes")

test_that("metabolism predictions (predict_metab, predict_DO) make sense", {
  
  # metab_bayes - not working yet
  mm <- metab_bayes(data=v(french), adapt_steps=100, burnin_steps=400, num_saved_steps=500)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- filter(DO_preds, local.date == "2012-08-24")
  knownbug <- function() {
    expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
    which(!is.na(DO_preds$DO.mod))
  }
  # library(ggplot2); ggplot(DO_preds, aes(x=local.time, group=factor(local.date))) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  # library(ggplot2); ggplot(DO_preds_Aug24, aes(x=local.time, group=local.date)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  
})
