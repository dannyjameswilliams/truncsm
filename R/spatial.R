#' Make Polygon Points
#'
#' @description Create a series of points in a polygon with a given resolution
#'
#' @param x vector of \code{x} vertices that form a closed cycle
#' @param y vector of \code{y} vertices that form a closed cycle
#' @param res resolution to set between vertices, i.e. number of points in each line segment
#'
#' @details \code{x} and \code{y} are vectors whose pairs give the vertices of the polygon object
#'
#' This is a simple function that takes the same initial two arguments as \code{\link[graphics]{polygon}}, but creates
#' points that go over the course of the polygon edges. To create a full shape, the final vertex of the pairs \code{x}
#' and \code{y} must match the first vertex.
#'
#' @return a matrix whose first column are the points in \code{x} and whose second column are the points in \code{y}
#' @export
#'
#' @examples
#' # make unit square
#'   polygon_points(c(-1,-1,1,1,-1),c(-1,1,1,-1,-1))
polygon_points = function(x, y, res = 20){
  out = cbind(seq(x[1],x[2], len=res),seq(y[1],y[2], len=res))
  for(j in 2:(length(x)-1)){
    out = rbind(out, cbind(seq(x[j],x[j+1], len=res),seq(y[j],y[j+1], len=res))[2:(res),])
  }
  return(out)
}


