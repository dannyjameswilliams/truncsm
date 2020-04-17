#' Euclidean distance on disk
#'
#' Boundary function which projects region onto 2D disk
#'
#' @param z matrix of euclidean coordinates
#' @param dV boundary in euclidean coordinates
#'
#' @details For data in Euclidean coordinates on a sphere, such as those from von-Mises Fisher distribution,
#' this function projects these points onto a 2D disk, and takes the Euclidean distance between the now 2D
#' points and the boundary.
#'
#' This is used for describing the boundary function \eqn{g}.
#'
#' @return List containing two elements
#' \item{\code{g}}{values of function \code{g}}
#' \item{\code{grad}}{matrix of derivatives of \code{g}}
#'
#' @export
g_disk = function(z, dV){
  grad = array(NA, dim(z))
  ztilde = array(NA, c(nrow(z), 2))
  for(i in 1:nrow(z)){
    diffs = dV[,1:2] - matrix(unlist(c(z[i,1:2])), nrow(dV), ncol(dV)-1, byrow=TRUE)
    ztilde[i,] = dV[which.min(sqrt(rowSums(diffs^2))), 1:2]
    P = diag(3) - z[i,] %*% t(z[i,])
    grad[i, ] = P %*% c((z[i,1:2] - ztilde[i,])/ sum(sqrt((z[i,1:2]-ztilde[i,])^2)), 0)
  }
  list("g"    = sqrt(rowSums((z[,1:2] - ztilde)^2)),
       "grad" = grad
  )
}

#' @keywords internal
hav_grad = function (x, y, r = 1) {
  .e1 <- (x[2] - y[2])/2
  .e2 <- (x[1] - y[1])/2
  .e3 <- cos(y[1])
  .e4 <- sin(.e1)
  .e5 <- cos(x[1])
  .e6 <- sin(.e2)
  .e7 <- .e4^2
  .e10 <- .e5 * .e3 * .e7 + .e6^2
  .e13 <- sqrt(1 - .e10) * sqrt(.e10)
  c(x1 = r * (cos(.e2) * .e6 - .e3 * .e7 * sin(x[1]))/.e13, x2 = r *
      cos(.e1) * .e5 * .e3 * .e4/.e13)
}

#' Haversine distance on sphere
#'
#' Calculate great circle distance between a set of points and a boundary, with gradient
#'
#' @param z matrix of euclidean coordinates
#' @param dV boundary in euclidean coordinates
#' @param r radius of sphere, default 1
#'
#' @return List containing two elements
#' \item{\code{g}}{values of function \code{g}}
#' \item{\code{grad}}{matrix of derivatives of \code{g}}
#'
#' @importFrom geosphere distHaversine
#' @export
g_hav = function(z, dV, r=1){
  x = euclid_to_sphere(z)
  dV = euclid_to_sphere(dV)
  dVpoints = (180/pi)*dV[,2:1] - cbind(rep(180, nrow(dV)), rep(90, nrow(dV)))

  grad = array(NA, dim(z))
  g = rep(NA, nrow(x))
  for(i in 1:nrow(x)){
    xpoint = (180/pi)*x[i,2:1] - c(180, 90)

    dists = distHaversine(xpoint, dVpoints, r=r)
    g[i] = min(dists)
    xtilde = dV[which.min(dists),]
    P = diag(3) - z[i,] %*% t(z[i,])
    grad[i,] = P %*% (g_derivx(z[i,]) %*% hav_grad(x[i,], xtilde, r=r))
  }
  list(g=g, grad=grad)
}


#' @keywords internal
get_J = function(z, dV, g, psi, r=1){

  # Calculate gradient in advance
  gvals = g(z, dV, r)
  g_grad = gvals$grad
  g_val = gvals$g

  # Now the gradients are defined in the get_J environment, the following function
  # will only have the given z, dV, and g in its environment, regardless whether they
  # are changed somwhere else. This also speeds up computation.
  J = function(mu){
    mu[1] = ifelse(mu[1] < 0, mu[1] + pi, mu[1])
    mu[1] = ifelse(mu[1] > pi, mu[1] - pi, mu[1])
    mu[2] = ifelse(mu[2] < 0, mu[2] + 2*pi, mu[2])
    mu[2] = ifelse(mu[2] > 2*pi, mu[2] - 2*pi, mu[2])

    mu = sphere_to_euclid(mu)

    Jout = rep(NA, nrow(z))
    for(i in 1:nrow(z)){
      psi_eval = psi(mu)
      grad = psi_eval$f; grad2 = psi_eval$grad
      P = diag(3) - z[i,] %*% t(z[i,])
      Jout[i] = g_val[i] * t(P %*% grad) %*% (P %*% grad) +
        2*g_val[i] * sum(diag(-2 * z[i,] %*% t(grad) + grad2 %*% P)) +
        2*(t(P %*% grad) %*% g_grad[i,])
    }
    return(0.5*mean(Jout))
  }
  return(J)
}

#' @keywords internal
g_default_sphere = function(z, dV){
  list("g" = rep(1, nrow(z)),
       "grad" = array(0, dim(dV))
       )
}

#' Truncated Spherical Score Matching
#'
#' @param x spherical coordinate data
#' @param dV boundary coordinates
#' @param family distribution, e.g. "vmf" for von-Mises Fisher
#' @param g boundary function
#' @param init optional; initial condition
#' @param options optional list; non-specified \code{family}
#' by including \code{psi} which returns first and second derivative of log pdf as \code{f} and \code{grad}
#'
#' @details if \code{g} is not included as an argument, then non-truncated score matching will be used. For
#' a sphere, use
#' \describe{
#' \item{\code{"Haversine"}}{haversine distance function \code{\link{g_hav}}}
#' \item{\code{"Euclidean"}}{euclidean projection distance function \code{\link{g_disk}}}
#' }
#' @return parameter estimates
#' @export
sphere_sm = function(x, dV, family="vmf", g, init=colMeans(euclid_to_sphere(x)), options = list()){
  psi = match_family_psi_sphere(family, options)
  g = match_g_sphere(g)
  obj_fn = get_J(x, dV, g, psi)
  est = optim(init, obj_fn, method="BFGS")
  return(est$par)
}

#' @keywords internal
match_g_sphere = function(g){
  types = c("Default", "Disk", "Euclidean", "Haversine")
  pick = grep(g, types, ignore.case = TRUE)
  if(pick == 1) return(g_default_sphere)
  if(pick == 2 | pick == 3) return(g_disk)
  if(pick == 4) return(g_hav)
  stop(paste0("g must be one of: ", types[1], ", ", types[3], ", ", types[4]))
}

#' @keywords internal
match_family_psi_sphere = function(name, options, psi=NULL){
  if(!is.null(options$psi)) psi = options$psi
  if(name=="vmf") psi = psi_vmf
  if(is.null(psi)) stop("either 'family' or 'psi' need to be specified, see ?sphere_sm")
  return(psi)
}

#' Derivatives of the Log Density of von-Mises Fisher distribution
#'
#' The first and second derivative of the log density of the von-Mises Fisher distribution
#'
#' @param mu mean direction to be estimated
#' @param x random variable x, taken one row at a time
#' @param k concentration parameter, default \code{k=10}
#'
#' @return list containing two elements
#' \item{\code{f}}{the evaluation of the first derivative}
#' \item{\code{grad}}{the evaluation of the second derivative}
#'
#' @export
psi_vmf = function(mu, x, k = 10) {
  return(list("f"=k*mu, "grad"=matrix(0, length(mu), length(mu))))
}



