[![Build Status](https://travis-ci.com/dannyjameswilliams/truncsm.svg?branch=master)](https://travis-ci.com/dannyjameswilliams/truncsm)  


# truncsm package

This package begins with an implementation of truncated score matching, developed by Song Liu et al, which you can find here:
[Estimating Density Models with Truncation Boundaries using Score Matching](https://arxiv.org/abs/1910.03834)

This work has been extended to score matching estimation on a manifold, with applications on a sphere. This is the R code accompaniment to the COMPASS mini project report. 

## Installation

**Currently the package is private, so cannot be installed in this way. See the vignettes for examples of the code being run**

This package can be installed through github, which you will need the `devtools` package for. To install `devtools`, [please see here](https://www.r-project.org/nosvn/pandoc/devtools.html). Once this is installed, you can install this package with the command in R:
```
devtools::install_github("dannyjameswilliams/truncsm")
```

## The Package

The aim of the package is to provide straightfoward access to estimating parameters in a truncated setting. The idea is to minimise the difference between the gradient of the log density for the model pdf and the data pdf, with a function that describes behaviour at the boundary of the objective function.

The functions `truncsm` and `sphere_sm` are the main feature of the package. Currently, you use the functions in the following way:

  1. Calculate `x`, the truncated dataset
  2. Form the boundary, i.e. a set of vertices that form a closed loop (so the first vertex = the second vertex)
  3. Transform this boundary into a high resolution set of points (you can use `polygon_points` in the case of the boundary being in the format of an R polygon)
  4. (Optional) set up functions for the first and second derivative of the log density, or the function that defines the boundary
     (Alternatively) specify the `family` argument
  5. (Optional) specify the scaling function. The default argument is no scaling function for `sphere_sm` and Euclidean distance function for `truncsm`. For `sphere_sm`, you may choose `Haversine` or `Projected Euclidean`
  5. Run `truncsm` or `sphere_sm` using these arguments, which will numerically optimise to find the parameter estimates.
  
 ## `sphere_sm` Examples
 Please see the html file (and Rmarkdown file) in `tutorials/` to see a few examples of estimating in the case of a 2D sphere, using the von-Mises Fisher distribution and the Kent distribution.

 ## `truncsm`: Multivariate Normal Example
 
 Here is some code which exemplifies the useage of this package to accurately estimate the mean of a multivariate Normal distribution centred on (2,2):
 
 ```
library(truncsm)
library(sp)

# sample mvn
n = 1000
real_mean = c(2,2)
real_cov = matrix(c(1,0,0,1),2,2)
x = cbind(rnorm(n, 2),rnorm(n, 2))

# get boundary and truncate points
xV = c(0,1,1,4,2,0)
yV = c(0.5,1,2,1,0,0.5)
dV = polygon_points(xV, yV)
truncated_points = point.in.polygon(x[,1],x[,2], xV, yV)
truncated_x = x[as.logical(truncated_points),]

# estimate mean
truncsm_est = truncsm(truncated_x, dV)

# plot
plot(x[!as.logical(truncated_points),1], x[!as.logical(truncated_points),2], 
     pch=20, xlab="x", ylab="y", xlim = range(x[,1])+c(0,1), 
     ylim=range(x[,2]) + c(0,1))
polygon(x=xV, y=yV, col=rgb(0,0,1,0.05), border=rgb(0,0,1,0.4), lwd=2)
points(truncated_x, col=rgb(0.2,0.2,0.7), pch=20)
points(truncsm_est[1], truncsm_est[2], bg="green", cex=2, pch=21)
points(2,2, bg="red", pch=21, cex=2)
 ```
