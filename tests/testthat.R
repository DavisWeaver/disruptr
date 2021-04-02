Sys.setenv("R_TESTS" = "")
library(testthat)
library(disruptr)

test_check("disruptr")
