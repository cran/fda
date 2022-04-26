## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# create a fine mesh of x-values over [0,2*pi]
x <- seq(0,2*pi,len=101)
# step size
delta <- 2*pi/100
# sin(x)
W <- sin(x)
# exponentiate the result
EW <- exp(W)
# compute the integral from 0 to 2*pi using the trapezoidal rule
hof2pi <- delta*(sum(EW) - EW[101]/2)
print(paste("h(2*pi) =",round(hof2pi,2)))

## ----fig.height = 7-----------------------------------------------------------
h <- delta*(cumsum(EW) - EW[101]/2)
par(mfcol=c(2,1))
plot(x, W, type="l")
plot(x, h, type="l")

## ----fig.height = 7-----------------------------------------------------------
h <- delta*(cumsum(EW) - EW[101]/2)/hof2pi
par(mfcol=c(2,1))
plot(x, W, type="l")
plot(x, h, type="l")

