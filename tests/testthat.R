Sys.setenv("R_TESTS" = "")
library(testthat)
test_check('streamMetabolizer', reporter="summary")
