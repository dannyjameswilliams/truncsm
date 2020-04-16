#' Translate Spherical Coordinates
#'
#' Switch between spherical coordinates \code{(theta, phi)} and euclidean coordinates \code{(z1, z2, z3)}
#'
#' @param x spherical input
#' @param z euclidean input
#'
#' @return vector or matrix of same dimension as the input
#' @export
#'
#' @examples
#' sph_x = c(pi/2, pi)
#' sph_z = sphere_to_euclid(sph_x)
#' sph_x_2 = euclid_to_sphere(sph_z)
sphere_to_euclid = function(x){
  xvector = is.vector(x)
  if(xvector) x = t(x)
  z1 = cos(x[,1])
  z2 = sin(x[,1])*cos(x[,2])
  z3 = sin(x[,1])*sin(x[,2])
  if(xvector) return(c(z1, z2, z3)) else return(cbind(z1, z2, z3))
}
#' @describeIn sphere_to_euclid
euclid_to_sphere = function(z){
  if(is.vector(z)) z = t(z)
  lon = atan2(z[,3], z[,2])
  r = sqrt(rowSums(z^2))
  colat = acos(z[,1]/r)
  if(is.vector(z)) {
    x=(c(colat, lon))
    x[1] = ifelse(x[1] < 0, x[1] + pi, x[1])
    x[2] = ifelse(x[2] < 0, x[2] + 2*pi, x[2])
    return(x)
  } else {
    x=(cbind(colat, lon))
    x[,1] = ifelse(x[,1] < 0, x[,1] + pi, x[,1])
    x[,2] = ifelse(x[,2] < 0, x[,2] + 2*pi, x[,2])
    return(x)
  }
}
