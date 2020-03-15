test_that("multiplication works", {
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
})
