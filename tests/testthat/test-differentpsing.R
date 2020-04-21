test_that("input different psi and g (NORMAL)", {
  psi_mvn = function(theta, x) {
    if(length(theta)>length(x)) Cov = matrix(theta[(length(x)+1):length(theta)],length(x),length(x)) else Cov = diag(length(theta))
    mu = theta[1:length(x)]
    Cov_inv = solve(Cov)
    return(list("f"=-Cov_inv %*% (x-mu), "grad"=-Cov_inv))
  }
  g_def = function(x, dV){
    xtilde = array(NA, dim(x))
    for(i in 1:nrow(x)) {
      diffs = dV-matrix(unlist(c(x[i,])),nrow(dV),ncol(dV),byrow=TRUE)
      xtilde[i,] = dV[which.min(sqrt(rowSums(diffs^2))),]
    }
    list("g"=sqrt(rowSums((x - xtilde)^2)),"grad"=(x - xtilde)/sqrt(rowSums((x-xtilde)^2)))
  }

  n = 1000
  real_mean = c(2,2)
  real_cov = matrix(c(1,0,0,1),2,2)
  x = cbind(rnorm(n, 2),rnorm(n, 2))
  x_bounds = c(0, 0, 1.5, 1.5, 0)
  y_bounds = c(0, 3, 3, 0, 0)
  dV = polygon_points(x_bounds, y_bounds, 50)
  library(sp)
  truncated_points = point.in.polygon(x[,1],x[,2], x_bounds, y_bounds)
  truncated_x = x[as.logical(truncated_points),]
  truncsm_est = truncsm(x=truncated_x, dV = dV, options=list(g=g_def, psi=psi_mvn))
  expect_true(1==1)
})

test_that("input different psi and g (SPHERICAL)", {
  psi_vmf2 = function(par, x, k=30){
    list(
      "f" = k*par,
      "grad" = matrix(0, nrow= length(par), ncol=length(par))
    )
  }
  g_def2 = function(x, dV, r=NULL){
    xtilde = array(NA, dim(x))
    for(i in 1:nrow(x)) {
      diffs = dV-matrix(unlist(c(x[i,])),nrow(dV),ncol(dV),byrow=TRUE)
      xtilde[i,] = dV[which.min(sqrt(rowSums(diffs^2))),]
    }
    list("g"=sqrt(rowSums((x - xtilde)^2)),"grad"=(x - xtilde)/sqrt(rowSums((x-xtilde)^2)))
  }

  n = 1000
  set.seed(1)
  real_mu = c(pi/2,pi)
  centre_euclid = sphere_to_euclid(real_mu)
  set.seed(1)
  x1 = rnorm(n, real_mu[1], 0.1/3)
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

  est = sphere_sm(z, z_dV, options = list(
    psi = psi_vmf2,
    g = g_def2
  ))
  expect_true(sqrt(sum((real_mu - est))^2) < 0.02)
})
