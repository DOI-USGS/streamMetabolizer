context("metab_Kmodel")

library(dplyr)
vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)

test_that("metab_Kmodel predictions (predict_metab, predict_DO) make sense", {
  
  # fit a first-round MLE and extract the K estimates
  expect_warning({mm1 <- metab_mle(data=vfrench, day_start=-1, day_end=23)}, "temperature out of range")
  K600_mm1 <- predict_metab(mm1) %>% select(local.date, K600, K600.lower, K600.upper)
  
  # smooth the K600s
  expect_warning({mm2 <- metab_Kmodel(data_daily=K600_mm1, method='mean', transforms=c(K600='log'))}, "no SE available")
  K600_mm2 <- predict_metab(mm2) %>% select(local.date, K600)
  
  # refit the MLE with fixed K
  expect_warning({mm3 <- metab_mle(data=vfrench, data_daily=K600_mm2, day_start=-1, day_end=23)}, "temperature out of range")
  predict_metab(mm3)
  
  # compare the two MLE results
  # library(ggplot2)
  # library(tidyr)
  # preds <- bind_rows(
  #   mutate(predict_metab(mm1), model='PRK'),
  #   mutate(predict_metab(mm3), model='PR')) %>%
  #   filter(!is.na(K600)) %>%
  #   gather(estimate, value, GPP, ER, K600)
  # ggplot(preds, aes(x=local.date, y=value, color=model)) + geom_point() + geom_line() + theme_bw() + facet_grid(estimate ~ ., scales = "free_y")
  
})