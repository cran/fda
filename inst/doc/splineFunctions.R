## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup fda package--------------------------------------------------------
library(fda)

## ----Two simple spline functions----------------------------------------------
par(mfrow=c(1,2))
# first graph: order 2
# set range, number of basis functions and order
rng    <- c(0,1) 
nbasis <- 3
norder <- 2
#  make the spline basis object
basis2 <- create.bspline.basis(rng, nbasis, norder)
#  define three coefficients
coefs2 <- matrix(c(1,2,-1), nbasis,1)
#  make the spline function object
splfd2 <- fd(coefs2, basis2)
#  plot the spline function with a vertical dashed 
#  line at the interior knot location
plot(splfd2, xlab="t", ylab="s(t)")
lines(c(0.5,0.5), c(-1,2), lwd=2, lty=2)
# second graph: order 3
# set range, number of basis functions and order
rng    <- c(0,1) 
nbasis <- 4
norder <- 3
#  make the spline basis object
basis3 <- create.bspline.basis(rng, nbasis, norder)
#  define four coefficients
coefs3 <- matrix(c(1,2,-1,3), nbasis,1)
#  make the spline function object
splfd3 <- fd(coefs3, basis3)
#  plot the spline function with a vertical dashed 
#  line at the interior knot location
plot(splfd3, xlab="t", ylab="s(t)", main="Order 3 spline function")
lines(c(0.5,0.5), c(0,3), lwd=2, lty=2)

