#  set up objects for examples

#  constant basis for estimating weight functions
cbasis = create.constant.basis(c(0,1))
#  monomial basis: {1,t}  for estimating weight functions
mbasis = create.monomial.basis(c(0,1),2)
#  quartic spline basis with 54 basis functions for
#    defining functions to be analyzed
xbasis = create.bspline.basis(c(0,1),24,5)
#  set up functional parameter objects for weight bases
cfdPar = fdPar(cbasis)
mfdPar = fdPar(mbasis)
#  sampling points over [0,1]
tvec = seq(0,1,len=101)

#  Example 1:  a single first order constant coefficient unforced equation
#     Dx = -4*x  for  x(t) = exp(-4t)

beta    = 4
xvec    = exp(-beta*tvec)
xfd     = smooth.basis(tvec, xvec, xbasis)$fd
xfdlist = list(xfd)
bwtlist = list(cfdPar)

result = pda.fd(xfdlist, bwtlist)

#  display weight coefficient for variable
bwtlistout = result$bwtlist
bwtfd      = bwtlistout[[1]]$fd
par(mfrow=c(1,1))
plot(bwtfd)
title("Weight coefficient for variable")
print(round(bwtfd$coefs,3))
#  display residual functions
reslist    = result$resfdlist
plot(reslist[[1]])
title("Residual function")

#  Example 2:  a single first order varying coefficient unforced equation
#     Dx(t) = -t*x(t) or x(t) = exp(-t^2/2)

bvec    = tvec
xvec    = exp(-tvec^2/2)
xfd     = smooth.basis(tvec, xvec, xbasis)$fd
xfdlist = list(xfd)
bwtlist = list(mfdPar)

result = pda.fd(xfdlist, bwtlist)

#  display weight coefficient for variable
bwtlistout = result$bwtlist
bwtfd      = bwtlistout[[1]]$fd
par(mfrow=c(1,1))
plot(bwtfd)
title("Weight coefficient for variable")
print(round(bwtfd$coefs,3))
#  display residual function
reslist    = result$resfdlist
plot(reslist[[1]])
title("Residual function")

#  Example 3:  a single second order constant coefficient unforced equation
#     Dx(t) = -(2*pi)^2*x(t) or x(t) = sin(2*pi*t)

xvec    = sin(2*pi*tvec)
xfd     = smooth.basis(tvec, xvec, xbasis)$fd
xfdlist = list(xfd)
bwtlist = list(cfdPar,cfdPar)

result = pda.fd(xfdlist, bwtlist)

#  display weight coefficients
bwtlistout = result$bwtlist
bwtfd1     = bwtlistout[[1]]$fd
bwtfd2     = bwtlistout[[2]]$fd
par(mfrow=c(2,1))
plot(bwtfd1)
title("Weight coefficient for variable")
plot(bwtfd2)
title("Weight coefficient for derivative of variable")
print(round(c(bwtfd1$coefs, bwtfd2$coefs),3))
print(bwtfd2$coefs)
#  display residual function
reslist    = result$resfdlist
par(mfrow=c(1,1))
plot(reslist[[1]])
title("Residual function")

#  Example 4:  two first order constant coefficient unforced equations
#     Dx1(t) = x2(t) and Dx2(t) = -x1(t)  
#   equivalent to  x1(t) = sin(2*pi*t)

N = 5
xmat1     = outer(sin(2*pi*tvec),rnorm(N))
xmat2     = outer(2*pi*cos(2*pi*tvec),rnorm(N))
xfd1      = smooth.basis(tvec, xvec1, xbasis)$fd
xfd2      = smooth.basis(tvec, xvec2, xbasis)$fd
xfdlist   = list(xfd1,xfd2)
bwtlist   = list(
                 list(
                      list(cfdPar),
                      list(cfdPar)
                     ),
                 list(
                      list(cfdPar),
                      list(cfdPar)
                     )
                )

result = pda.fd(xfdlist, bwtlist)

#  display weight coefficients
bwtlistout = result$bwtlist
bwtfd11    = bwtlistout[[1]][[1]][[1]]$fd
bwtfd21    = bwtlistout[[2]][[1]][[1]]$fd
bwtfd12    = bwtlistout[[1]][[2]][[1]]$fd
bwtfd22    = bwtlistout[[2]][[2]][[1]]$fd
par(mfrow=c(2,2))
plot(bwtfd11)
title("Weight for variable 1 in equation 1")
plot(bwtfd21)
title("Weight for variable 2 in equation 1")
plot(bwtfd12)
title("Weight for variable 1 in equation 2")
plot(bwtfd22)
title("Weight for variable 2 in equation 2")
print(round(bwtfd11$coefs,3))
print(round(bwtfd21$coefs,3))
print(round(bwtfd12$coefs,3))
print(round(bwtfd22$coefs,3))
#  display residual functions
reslist = result$resfdlist
par(mfrow=c(2,1))
plot(reslist[[1]])
title("Residual function for variable 1")
plot(reslist[[2]])
title("Residual function for variable 2")

#  Example 5:  a single first order constant coefficient equation
#     Dx = -4*x  for  x(t) = exp(-4t) forced by u(t) = 2

beta    = 4
alpha   = 2
xvec0   = exp(-beta*tvec)
intv    = (exp(beta*tvec) - 1)/beta
xvec    = xvec0*(1 + alpha*intv)
xfd     = smooth.basis(tvec, xvec, xbasis)$fd
xfdlist = list(xfd)
bwtlist = list(cfdPar)
awtlist = list(cfdPar)
ufdlist = list(fd(1,cbasis))

result = pda.fd(xfdlist, bwtlist, awtlist, ufdlist)

#  display weight coefficients
bwtlistout = result$bwtlist
bwtfd      = bwtlistout[[1]]$fd
awtlistout = result$awtlist
awtfd      = awtlistout[[1]]$fd
par(mfrow=c(2,1))
plot(bwtfd)
title("Weight for variable")
plot(awtfd)
title("Weight for forcing function")
#  display residual function
reslist = result$resfdlist
par(mfrow=c(1,1))
plot(reslist[[1]], ylab="residual")
title("Residual function")

#  Example 6:  two first order constant coefficient equations
#     Dx = -4*x    for  x(t) = exp(-4t)     forced by u(t) =  2
#     Dx = -4*t*x  for  x(t) = exp(-4t^2/2) forced by u(t) = -1

beta    = 4
xvec10  = exp(-beta*tvec)
alpha1  = 2
alpha2  = -1
xvec1   = xvec10 + alpha1*(1-xvec10)/beta
xvec20  = exp(-beta*tvec^2/2)
vvec    = exp(beta*tvec^2/2);
intv    = 0.01*(cumsum(vvec) - 0.5*vvec)
xvec2   = xvec20*(1 + alpha2*intv)
xfd1    = smooth.basis(tvec, xvec1, xbasis)$fd
xfd2    = smooth.basis(tvec, xvec2, xbasis)$fd
xfdlist = list(xfd1, xfd2)
bwtlist    = list(
                 list(
                      list(cfdPar),
                      list(cfdPar)
                     ),
                 list(
                      list(cfdPar),
                      list(mfdPar)
                     )
                )
awtlist = list(list(cfdPar), list(cfdPar))
ufdlist = list(list(fd(1,cbasis)), list(fd(1,cbasis)))

result = pda.fd(xfdlist, bwtlist, awtlist, ufdlist)

# display weight functions for variables
bwtlistout = result$bwtlist
bwtfd11    = bwtlistout[[1]][[1]][[1]]$fd
bwtfd21    = bwtlistout[[2]][[1]][[1]]$fd
bwtfd12    = bwtlistout[[1]][[2]][[1]]$fd
bwtfd22    = bwtlistout[[2]][[2]][[1]]$fd
par(mfrow=c(2,2))
plot(bwtfd11)
title("weight on variable 1 in equation 1")
plot(bwtfd21)
title("weight on variable 2 in equation 1")
plot(bwtfd12)
title("weight on variable 1 in equation 2")
plot(bwtfd22)
title("weight on variable 2 in equation 2")
print(round(bwtfd11$coefs,3))
print(round(bwtfd21$coefs,3))
print(round(bwtfd12$coefs,3))
print(round(bwtfd22$coefs,3))
#  display weight functions for forcing functions
awtlistout = result$awtlist
awtfd1     = awtlistout[[1]][[1]]$fd
awtfd2     = awtlistout[[2]][[1]]$fd
par(mfrow=c(2,1))
plot(awtfd1)
title("weight on forcing function in equation 1")
plot(awtfd2)
title("weight on forcing function in equation 2")
#  display residual functions
reslist    = result$resfdlist
par(mfrow=c(2,1))
plot(reslist[[1]])
title("residual function for equation 1")
plot(reslist[[2]])
title("residual function for equation 2")



