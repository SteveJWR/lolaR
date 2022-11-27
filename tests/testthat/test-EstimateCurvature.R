test_that("Verify Vectorized Midpoint Computation", {
  x = seq(-0.001,0.001,length.out = 10001)
  expect_vector(MidDist(x,1,1,1), ptype = double(), size = 10001)
})

test_that("Verify Curvature Estimation Works", {
  expect_equal(EstimateKappa(1,1,1,sqrt(3/4)), 0)
})


