library(testthat)
library(onlineBcp)

x <- c(rnorm(10, 0, 1), rnorm(10, 5, 1))
bcp <- online_cp(x)
summary(bcp)


test_check("onlineBcp")
