test_that("polygons work", {
  polygon_points(c(-1,-1,1,1,-1),c(-1,1,1,-1,-1))
  polygon_points(c(-1,-1,1,1,2,2,3,3,4,4,5,5),c(-1,1,1,-1,c(2,2,3,3,4,4,5,5)*-1))
})
