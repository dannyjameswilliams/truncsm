library(sp)

test_that("mvn simulation", {
  mvn_sim = function(){
    n = 1000
    real_mean = c(2,2)
    real_cov = matrix(c(1,0,0,1),2,2)
    x = cbind(rnorm(n, 2),rnorm(n, 2))
    x_bounds = c(0, 0, 1.5, 1.5, 0)
    y_bounds = c(0, 3, 3, 0, 0)
    dV = polygon_points(x_bounds, y_bounds, 50)
    truncated_points = point.in.polygon(x[,1],x[,2], x_bounds, y_bounds)
    truncated_x = x[as.logical(truncated_points),]
    truncsm_est = truncsm(x=truncated_x, dV = dV, family=mvn(Cov = 1))
    return(norm(real_mean - truncsm_est, "2"))
  }
  out = rep(NA, 100)
  for(i in 1:100) {
    set.seed(i)
    out[i] = mvn_sim()
  }
  expect_true(mean(out) < 0.2)
})
