
require(testthat)
require(rtnc)

fn <- function(x){
  x <- as.double(x)
  f = x[1]^2 +  abs(x[2])^2
  g = c(0,0)
  g[1] = 2 * x[1]
  g[2] = 3 * abs(x[2])^2
  g[2] = ifelse(x[2] < 0, -g[2], g[2])
  return(list(f, g))
}


test_that("The replication of standard example works",{
  res = tnc(c(1,2), fn, maxCGit=2, maxnfeval=20, low=c(NA, 1))
  expect_that(res$x, equals(c(0,1)))
  expect_that(res$f, equals(1))
  expect_that(res$g, equals(c(0,3)))
})
