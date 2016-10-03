context('distribs')
# this file is really just to help me wrap my head around parameter
# distributions and parameter scaling.

library(ggplot2)

# variable description convention in notes: 
# stan.and.SM.arg = symbol = r.arg

test_that("scaling from standard normal to lognormal", {
  skip('tests are for developer understanding, not code performance')
  
  # location = mu = meanlog; scale = sigma = sdlog
  location <- 0
  scale <- 1
  n <- 10000
  
  # desired
  desired <- rlnorm(n, meanlog=location, sdlog=scale)
  expect_equal(mean(log(desired)), location, tol=0.5)
  expect_equal(sd(log(desired)), scale, tol=0.2)
  expect_gt(min(desired), 0)
  expect_equal(min(desired), 0, tol=0.1)
  ggplot(data_frame(desired), aes(x=log(desired))) + geom_histogram()
  ggplot(data_frame(desired), aes(x=desired)) + geom_histogram() + xlim(-1, 5)
  
  # scaled
  scaled <- rnorm(n, mean=0, sd=1)
  ggplot(data_frame(scaled), aes(x=scaled)) + geom_histogram()
  final <- exp(location)*(exp(scaled)^scale)
  ggplot(data_frame(final), aes(x=log(final))) + geom_histogram()
  ggplot(data_frame(final), aes(x=final)) + geom_histogram() + xlim(-5, 100)
  
  # comparison
  comp <- bind_rows(data_frame(value=desired, type='desired'), data_frame(value=scaled, type='scaled'), data_frame(value=final, type='final'))
  ggplot(comp) + geom_histogram(aes(x=log(value), fill=type), alpha=0.5, binwidth=1, position='identity')
  ggplot(filter(comp, type!='scaled')) + geom_histogram(aes(x=value, fill=type), alpha=0.5, position='identity') + xlim(-5, 100)
})
