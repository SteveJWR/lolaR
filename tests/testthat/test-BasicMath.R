test_that("secant works", {
  expect_equal(sec(0), 1)
})

test_that("expit works", {
  expect_equal(expit(0), 1/2)
})

test_that("SoftThreshold works", {
  x <- seq(-10,10,length.out = 100)
  expect_vector(SoftThreshold(x,10))
})


