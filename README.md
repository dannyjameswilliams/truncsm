# truncsm package

This package is an implementation of truncated score matching, developed by Song Liu et al, which you can find here:
[Estimating Density Models with Complex Truncation Boundaries](https://arxiv.org/abs/1910.03834)

**This package is a work in progress.**

## Installation

This package can be installed through github, which you will need the `devtools` package for. To install `devtools`, [please see here](https://www.r-project.org/nosvn/pandoc/devtools.html). Once this is installed, you can install this package with the command in R:
```
devtools::install_github("dannyjameswilliams/truncsm")
```

## The Package

The aim of the package is to provide straightfoward access to estimating parameters in a truncated setting. The idea is to minimise
the difference between the gradient of the log density for the model pdf and the data pdf, with a function that describes
behaviour at the boundary of the objective function.

The function `truncsm` is the main feature of the package. Currently, you use the function in the following way:

  1. Calculate `x`, the truncated dataset
  2. Form the boundary, i.e. a set of vertices that form a closed loop (so the first vertex = the second vertex)
  3. Use the function `polygon_points` to transform this boundary into a sequence of points
  4. (Optional) set up functions for the first and second derivative of the log density, or the function that defines the boundary
     (Alternatively) specify the `family` argument, however this is currently only available for the multivariate normal distribution.
  5. Run `truncsm` using these arguments, which will numerically optimise to find the parameter estimates.
  

 ## Multivariate Normal Example
 
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
