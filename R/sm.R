#' Default Boundary Function
#'
#' Default boundary function \eqn{g_0(x)} as described in "Estimating Density Models with Complex Truncation Boundaries"
#' by Song Liu et al
#'
#' @param x single point in truncated region
#' @param dV boundary in points
#'
#' @return list with two elements
#' \item{\code{g}}{values of \code{g}}
#' \item{\code{grad}}{matrix of derivatives of \code{g}}
#'
#' @export
g_def = function(x, dV){
  xtilde = array(NA, dim(x))
  for(i in 1:nrow(x)) {
    diffs = dV-matrix(unlist(c(x[i,])),nrow(dV),ncol(dV),byrow=TRUE)
    xtilde[i,] = dV[which.min(sqrt(rowSums(diffs^2))),]
  }
  list("g"=sqrt(rowSums((x - xtilde)^2)),"grad"=(x - xtilde)/sqrt(rowSums((x-xtilde)^2)))
}

#' Objective Function
#' @keywords internal
J_ = function(theta, x, dV, psi, g_val, g_grad){
  dpsi_eval = array(NA, c(nrow(x),ncol(dV),ncol(dV)))
  psi_eval = array(NA, c(nrow(x), ncol(dV)))
  for(i in 1:nrow(x)) {
    psi_out = psi(theta, x[i,])
    dpsi_eval[i,,] = psi_out$grad
    psi_eval[i,] = psi_out$f
  }
  term1 = mean(rowSums(psi_eval^2 * g_val))
  term2 = 2*mean(sum(apply(dpsi_eval, 1, function(y)sum(diag(y)))*g_val))
  term3 = 2*mean(rowSums(psi_eval*g_grad))
  term1 + term2 + term3
}

#' Truncated Score Matching
#'
#' For truncated data within some boundary, estimate parameters that do not necessarily rely on the boundary.
#'
#' @param x truncated random observation
#' @param dV boundary where \code{x} was truncated (in points)
#' @param family distribution family to look in, currently supports multivariate normal as "mvn"
#' @param init initial conditions to be passed to \code{\link[stats]{optim}}
#' @param options extra options that can be specified, see details
#'
#' @details
#' For details of the procedure, see "Estimating Density Models with Complex Truncation Boundaries" by Song Liu et al.
#'
#' The variable \code{x} takes values within a certain truncated region, where the boundary of this region is defined by
#' \code{dV}. It is possible through score matching to estimate parameters that are not constrained
#' to the boundary by minimising the difference in score functions of the model and the data - the gradient of the log pdf, i.e.
#' \deqn{
#' \psi = \nabla log p(x; \theta).
#' }
#'
#' If \code{family} is not supplied, then \code{psi} needs to be supplied as an element of \code{options}.
#' \code{psi} will be function will calculates the derivatives of the log pdf, and needs to output a list
#' containing two elements: the first derivative as \code{f} and second derivative as \code{grad}.
#'
#' All arguments to \code{options} are:
#' \describe{
#' \item{\code{psi}}{function outputting first and second derivatives of log pdf as \code{f} and \code{grad}, taking two inputs only: \code{theta} and \code{x}; the parameter to estimate over and the truncated observations.}
#' \item{\code{g}}{function outputting function controlling behaviour at the boundary, outputting value and first derivative of the function as \code{g} and \code{grad}, taking two inputs only: the truncated observations, \code{x} and the boundary points \code{dV}.}
#' }
#' See \code{\link{psi_mvn}} and \code{\link{g_def}} for examples.
#'
#' @importFrom stats optim
#' @return parameter estimates given by score matching
#' @export
truncsm = function(x, dV, family="mvn", init=rep(0.1,ncol(x)),
                   options = list()){
  psi = match_family_psi(family, options)
  g = match_g(options)
  g_eval = g(x, dV)
  out = optim(init, J_, x=x, dV=dV, psi=psi, g_val=g_eval$g, g_grad = g_eval$grad,
              method = "BFGS")
  return(out$par)
}

#' @keywords internal
match_g = function(options){
  if(is.null(options$g)) g = g_def else g=options$g
  return(g)
}

#' @keywords internal
match_family_psi = function(name, options, psi=NULL){
  if(!is.null(options$psi)) {
    return(options$psi)
  }
  if(name=="mvn") return(psi_mvn)
  stop("either 'family' or 'psi' need to be specified, see ?truncsm")
}

#' Derivatives of the Log Density of Multivariate Normal Distribution
#'
#' The first and second derivative of the log density of the Multivariate Normal Distribution
#'
#' @param theta parameters to be estimated, can include elements of the covariance matrix
#' @param x random variable x, taken one row at a time
#'
#' @return list containing two elements
#' \item{\code{f}}{the evaluation of the first derivative}
#' \item{\code{grad}}{the evaluation of the second derivative}
#'
#' @export
psi_mvn = function(theta, x) {
  if(length(theta)>length(x)) Cov = matrix(theta[(length(x)+1):length(theta)],length(x),length(x)) else Cov = diag(length(theta))
  mu = theta[1:length(x)]
  Cov_inv = solve(Cov)
  return(list("f"=-Cov_inv %*% (x-mu), "grad"=-Cov_inv))
}
