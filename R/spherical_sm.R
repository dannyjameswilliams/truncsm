#' Euclidean distance on disk
#'
#' Boundary function which projects region onto 2D disk
#'
#' @param z matrix of euclidean coordinates
#' @param dV boundary in euclidean coordinates
#' @param r radius of sphere
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
g_disk = function(z, dV, r = NULL){
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
#' @details This is one possible boundary function \code{g} that can be used for truncated score matching on a sphere,
#' used by the \code{"Haversine"} argument to \code{\link{sphere_sm}}.
#'
#' This function first finds the distance between the truncated dataset \code{z} and the boundary \code{dV}, using the
#' haversine distance function (from the \code{geosphere} package):
#' \deqn{
#'  d = 2 r arcsin( \sqrt{ sin^2 ( (\theta_2 - \theta_1) / 2 ) + cos (\theta_1) cos(\theta_2) sin^2( (\phi_2 - \phi_1)/2 ) } ),
#' }
#' where \eqn{(\phi_1, \theta_1)} is the longitude and latitude of point 1, and \eqn{(\phi_2, \theta_2)}
#' is the longitude and latitude of point 2, both in radians, and \eqn{r} is the radius, default \code{r=1}.
#' Once the distances are found, the smallest distances between \code{x} and the boundary as well as the first derivative
#' of the haversine function are saved and output back to the score matching function.
#'
#' Note that \code{z} and \code{dV} are Euclidean coordinates, and this function converts them to spherical coordinates
#' before taking the distance (as required by the Haversine). This is for consistency.
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
  # are changed somewhere else. This also speeds up computation.
  J = function(mu){
    if(length(mu) > ncol(z)) dim(mu) = c(ncol(z)-1, length(mu)/(ncol(z)-1))
    mu = spherical_bounds(mu)
    mu = par_to_euclid(mu)

    Jout = rep(NA, nrow(z))
    for(i in 1:nrow(z)){
      z0 = unlist(z[i,])
      psi_eval = psi(mu, z0)
      grad = psi_eval$f; grad2 = psi_eval$grad
      P = diag(3) - z0 %*% t(z0)
      Jout[i] = g_val[i] * t(P %*% grad) %*% (P %*% grad) +
        2*g_val[i] * sum(diag(-2 * z0 %*% t(grad) + grad2 %*% P)) +
        2*(t(P %*% grad) %*% g_grad[i,])
    }
    return(0.5*mean(Jout))
  }
  return(J)
}

#' @keywords internal
g_default_sphere = function(z, dV, r=NULL){
  list("g" = rep(1, nrow(z)),
       "grad" = array(0, dim(z))
       )
}

#' Truncated Spherical Score Matching
#'
#' Numerical minimisation of objective function corresponding to score matching on a sphere with a boundary
#' defined on a sphere
#'
#' @param x (truncated) euclidean coordinate data; 3D
#' @param dV boundary euclidean coordinates; 3D
#' @param family distribution, e.g. "vmf" for von-Mises Fisher
#' @param g boundary function
#' @param init optional; initial condition
#' @param options optional list; non-specified \code{family}
#' by including \code{psi} which returns first and second derivative of log pdf as \code{f} and \code{grad}
#'
#' @details
#' The variable \code{x} can take values within a certain truncated region, where the boundary of this region is defined by
#' \code{dV}. It is possible through score matching to estimate parameters that are not constrained
#' to the boundary by minimising the difference in score functions of the model and the data - the gradient of the log pdf, i.e.
#' \deqn{
#' \psi = \nabla log p(x; \theta).
#' }
#'
#' This function is mostly a wrapper function; it goes through multiple functions to eventually optimise
#' the objective function. If \code{g} is not included as an argument, then non-truncated score matching will be used. For
#' a sphere.
#'
#' Use the following arguments for \code{g}:
#' \describe{
#' \item{\code{"Haversine"}}{haversine distance function \code{\link{g_hav}}}
#' \item{\code{"Euclidean"}}{euclidean projection distance function \code{\link{g_disk}}}
#' }
#' The argument \code{g} can be overwritten by including a function as an element of the \code{options} list. This
#' needs to return the values of the function for each point in the truncated dataset as \code{g}, and its first derivative
#' as \code{grad}, both in a list. See \code{\link{g_hav}} and \code{\link{g_disk}} for examples.
#'
#' \code{psi} can also be overwritten by including it as an element of the \code{options} list. This needs to
#' return the first and second derivative of the log pdf as \code{f} and \code{grad} respectively in a list. See \code{\link{psi_vmf}}
#' or \code{\link{psi_kent}} for examples.
#' @return parameter estimates
#' @export
sphere_sm = function(x, dV, family="vmf", g="Default", init=colMeans(euclid_to_sphere(x)), options = list()){
  x = as.matrix(x)
  psi = match_family_psi_sphere(family, options)
  g_fn = match_g_sphere(g, options)
  obj_fn = get_J(x, dV, g_fn, psi)
  est = optim(init, obj_fn, method="BFGS")
  return(est$par)
}

#' @keywords internal
match_g_sphere = function(g, options){
  if(!is.null(options$g)) {
    return(options$g)
  }
  types = c("Default", "Disk", "Euclidean", "Haversine")
  pick = grep(g, types, ignore.case = TRUE)
  if(pick == 1) return(g_default_sphere)
  if(pick == 2 | pick == 3) return(g_disk)
  if(pick == 4) return(g_hav)
  stop(paste0("g must be one of: ", types[1], ", ", types[3], ", ", types[4]))
}

#' @keywords internal
match_family_psi_sphere = function(name, options){
  if(!is.null(options$psi)) {
    psi = options$psi
    return(psi)
  }
  types = c("vmf", "von Mises Fisher", "von-Mises Fisher", "kent")
  pick = grep(name, types, ignore.case = TRUE)
  if(pick==1 | pick==2 | pick==3) return(psi_vmf)
  if(pick==4) return(psi_kent)
  if(is.null(psi)) stop("either 'family' or 'psi' need to be specified, see ?sphere_sm")
}

#' Derivatives of the Log Density of spherical distributions
#'
#' The first and second derivative of the log density of the von-Mises Fisher distribution and Kent distribution
#'
#' @param mu mean direction to be estimated for the von-Mises Fisher distribution
#' @param par matrix of parameters to estimate for the Kent distribution
#' @param x random variable x, taken one row at a time
#' @param k concentration parameter, default \code{k=10}
#' @param beta ovalness parameter, default \code{beta=1}
#'
#' @details \code{psi_vmf} corresponds to the von-Mises Fisher distribution, where \code{psi_kent} corresponds to the
#' Kent distribution. The parameter to be estimated for the von-Mises Fisher distribution is the mean direction \code{mu},
#' whereas for the Kent distribution, the parameters are \code{mu}, and the major and minor axes. The parameters enter
#' \code{psi_kent} in a matrix with columns corresponding to these parameters respectively.
#'
#' @return list containing two elements
#' \item{\code{f}}{the evaluation of the first derivative}
#' \item{\code{grad}}{the evaluation of the second derivative}
#'
#' @name spherical_dist
NULL

#' @rdname spherical_dist
#' @export
psi_vmf = function(mu, x, k = 10) {
  return(list("f"=k*mu, "grad"=matrix(0, length(mu), length(mu))))
}

#' @rdname spherical_dist
#' @export
psi_kent = function(par, x, k=10, beta=1){
  gamma1 = par[,1]
  gamma2 = par[,2]
  gamma3 = par[,3]
  list(
    "f" = k*gamma1 + 2*beta*(gamma2 %*% t(gamma2 %*% x) - gamma3 %*% t(gamma3 %*% x)),
    "grad" = 2*beta*(gamma2 %*% t(gamma2) - gamma3 %*% t(gamma3))
  )
}



