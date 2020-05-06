#' @keywords internal
par_to_euclid = function(x){
  times = if(length(x) >= 6) 6 else 2
  xout = c()
  for(i in seq(1, times, by = 2)){
    xout = append(xout, cos(x[i]))
    xout = append(xout, sin(x[i])*cos(x[i+1]))
    xout = append(xout, sin(x[i])*sin(x[i+1]))
  }
  if(length(x) > times) xout = append(xout,
                                             x[(times+1):length(x)])
  return(xout)
}


#' Convert Spherical/Euclidean Coordinates
#'
#' @param x either spherical or Euclidean input
#' @name convert
NULL

#' @rdname convert
#' @export
sphere_to_euclid = function(x){
  xvector = is.vector(x)
  if(xvector) x = t(x)
  r = if(ncol(x) == 3) x[1,3] else 1
  z1 = r*cos(x[,1])
  z2 = r*sin(x[,1])*cos(x[,2])
  z3 = r*sin(x[,1])*sin(x[,2])
  if(xvector){
    return(c(z1, z2, z3))
  } else {
    return(cbind(z1, z2, z3))
  }
}

#' @rdname convert
#' @export
euclid_to_sphere = function(x){
  if(is.vector(x)) x = t(x)
  lon = atan2(x[,3], x[,2])
  r = sqrt(rowSums(x^2))
  colat = acos(x[,1]/r)
  if(is.vector(x)) {
    xout=(c(colat, lon))
    xout[1] = ifelse(xout[1] < 0, xout[1] + pi, xout[1])
    xout[2] = ifelse(xout[2] < 0, xout[2] + 2*pi, xout[2])
    return(xout)
  } else {
    xout=(cbind(colat, lon))
    xout[,1] = ifelse(xout[,1] < 0, xout[,1] + pi, xout[,1])
    xout[,2] = ifelse(xout[,2] < 0, xout[,2] + 2*pi, xout[,2])
    return(xout)
  }
}

#' @keywords internal
g_derivx = function(z){
  sqsum = z[1]^2 + z[2]^2 + z[3]^2
  sqsumxy = z[2]^2 + z[3]^2
  matrix(c(
    0, -sqrt(sqsumxy)/sqsum,
    -z[3]/sqsumxy, (z[2]*z[1])/(sqrt(sqsumxy)*sqsum),
    z[2]/sqsumxy, (z[3]*z[1])/(sqrt(sqsumxy)*sqsum)
  ),3,2,byrow=TRUE)[,c(2,1)]
}

#' @keywords internal
par_bounds = function(x){
  times = if(length(x) == 6) 6 else 2
  xout = c()
  for(i in seq(1, times, by = 2)){
    x[i] = x[i] %% (pi)
    x[i+1] = x[i+1] %% (2*pi)
  }
  return(x)
}

#' Convert Degrees to radian
#'
#' Take array of Colatitude, Longitude in degrees and convert to radians
#' @param x input in degrees
#'
#' @return \code{x} but in radians instead
radian = function(x){
  x*pi/180
}





