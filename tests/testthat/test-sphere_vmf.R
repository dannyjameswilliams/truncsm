test_that("estimate VMF mu on hemisphere", {
  n = 1000
  real_mu = c(pi/2,pi)
  centre_euclid = sphere_to_euclid(real_mu)
  set.seed(1)
  x1 = rnorm(n, real_mu[1], 0.1)
  x1 = ifelse(x1 > pi, x1-pi, x1)
  x1 = ifelse(x1 < 0, x1+pi, x1)
  x2 = rnorm(n, real_mu[2], 0.1)
  x2 = ifelse(x2 > 2*pi, x2-2*pi, x2)
  x2 = ifelse(x2 < 0, x2+2*pi, x2)
  x = cbind(x1, x2)
  z = sphere_to_euclid(x)

  boundary = pi/2
  dV = cbind(rep(pi/2, 100), seq(0, 2*pi, len=100))
  z_dV = sphere_to_euclid(dV)

  above_boundary = x[,1] > pi/2
  truncated_x = x[above_boundary,]
  truncated_z = z[above_boundary,]

  est = sphere_sm(z, z_dV, g="Haversine")
  expect_true(sqrt(sum((real_mu - est))^2) < 0.04)
})
