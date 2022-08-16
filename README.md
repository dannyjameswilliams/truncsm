# truncsm package

This package begins with my implementation of truncated score matching, or *TruncSM*, which is a method developed by Song Liu, Takafumi Kanamori and myself, which you can find here:
[Estimating Density Models with Truncation Boundaries using Score Matching](https://arxiv.org/abs/1910.03834)<sup>1</sup>

This work has been extended to score matching estimation on a manifold, with applications on a sphere, which is also implemented in this package. This work was developed by myself and Song Liu, and can be found here: [Score Matching for Truncated Density Estimation on a Manifold
](https://arxiv.org/abs/2206.14668)<sup>2</sup>

## Installation

This package can be installed through github, which you will need the `devtools` package for. To install `devtools`, [please see here](https://www.r-project.org/nosvn/pandoc/devtools.html). Once this is installed, you can install this package with the command in R:
```
devtools::install_github("dannyjameswilliams/truncsm")
```

## Overview 

The aim of the R package is to provide straightfoward access to estimating parameters in a truncated setting. The idea is to minimise the difference between the gradient of the log density for the model pdf and the data pdf, with a function that describes behaviour at the boundary of the objective function.

Please see the attached papers for details on the implementation. In its current state, you can use `truncsm` and `sphere_sm` using mostly default arguments.

The functions `truncsm` and `sphere_sm` are the main feature of the package. Currently, you use the functions in the following way:

  1. Given a truncated dataset, `X`, whose truncation is caused by a physical boundary, such as a series of lines (e.g. a country's borders) or otherwise
  2. Given a set of boundary points, `dV`, which make up a high resolution point set of the physical boundary that has truncated the dataset `X` (e.g. longitude and latitude of the country's borders)
  3. (Optional) set up functions for the first and second derivative of the log density, or the function that defines the boundary
     (Alternatively) specify the `family` argument
  4. (Optional) specify the scaling function. The default argument is no scaling function for `sphere_sm` and Euclidean distance function for `truncsm`. For `sphere_sm`, you may choose `Haversine` or `Projected Euclidean`.
  5. Run `truncsm` or `sphere_sm` using these arguments, which will numerically optimise to find the parameter estimates.

## Manifold Score Matching Examples

### Main Tutorial 
Please see the Rmarkdown file in `tutorials/` to see a few examples of estimating in the case of a 2D sphere, using the von-Mises Fisher distribution and the Kent distribution. To run this file, make sure this package is installed (see above), then open `tutorials/simulate.Rmd` in Rstudio. Provided you have the necessary installations for rendering Rmarkdown content, you should be able to knit the document to a HTML or PDF, and view the results there.

### Quick Example
A quick example is given in the code snippet below. First, simulate some data from a von-Mises Fisher distribution.
 
```r
library(truncsm)
n = 1000              
centre = c(pi/2, pi) 
centre_euclid = sphere_to_euclid(centre)   
set.seed(4)
z = Rfast::rvmf(n, centre_euclid, k = 6)
x = euclid_to_sphere(z)

dV = cbind(rep(pi/2+0.2, 200), seq(0, 2*pi, len=200))
outside_boundary = x[,1] > (pi/2+0.2)
truncated_x = x[outside_boundary,]

```
then estimate using `sphere_sm` (assuming the concentration parameter is known at `k=6`)
```r
est = sphere_sm(truncated_x, dV, g = "Haversine", family = vmf(k=6))
```
then you can plot 
```r
library(ggplot2)
ggplot() + geom_point(aes(x=x[,1], y=x[,2]), col="grey")+ geom_point(aes(x=truncated_x[,1], y=truncated_x[,2])) +
  geom_point(aes(est_hav$mu[1], est_hav$mu[2]), col="red", size=4) + theme_minimal()
```
and you should expect:

<img src="https://user-images.githubusercontent.com/56155783/173840452-9190b539-0b49-4d85-b072-36b651123f28.png" alt="Example estimate" width="500"/>

 ## `truncsm`: Multivariate Normal Example

 Here is some code which exemplifies the useage of this package to accurately estimate the mean of a multivariate Normal distribution centred on (2,2):
 
 ```r
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
 
 ## References
 
 [1] Liu, S., Kanamori, T., and Williams, D. J. Estimating density models with truncation boundaries using score matching. arXiv preprint, arXiv:1910.03834, 2022.
 
 [2] Williams, D. J. and Liu, S., Score Matching for Truncated Density Estimation on a Manifold. arXiv preprint, arXiv:2206.14668, 2022.
