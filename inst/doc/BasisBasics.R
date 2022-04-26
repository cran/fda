## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(fda)

## ----seven fourier basis functions, echo=FALSE--------------------------------
# set up a fourier basis object having 7 basis functions
nbasis <- 7
period <- 1
rangeval <- c(0,1)
fourier.basis <- create.fourier.basis(rangeval,nbasis, period)
# set up a plot of all of the seven basis functions
coefs <- diag(rep(1,nbasis))
coefs[1,1] <- 1
fourier.fd <- fd(coefs, fourier.basis)
plot(fourier.fd, xlab="t", ylab="phi(t)")

## ----one randomm fourier function---------------------------------------------
coefs <- matrix(rnorm(7),7,1)
onefourier.fd <- fd(coefs, fourier.basis)
plot(onefourier.fd, xlab="t", ylab="f(t)")

## ----spline basis functions---------------------------------------------------
nbasis <- 8
rangeval <- c(0,10)
spline.basis <- create.bspline.basis(rangeval, nbasis)
coefs <- diag(rep(1,8))
spline.fd <- fd(coefs, spline.basis)
plot(spline.fd, xlab = "t", ylab = "phi(t)")

## ----one random spline function-----------------------------------------------
coefs <- matrix(rnorm(24),8,3)
threespline.fd <- fd(coefs, spline.basis)
plot(threespline.fd, xlab="t", ylab="f(t)")

