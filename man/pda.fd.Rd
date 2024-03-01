\name{pda.fd}
\alias{pda.fd}
\title{
  Principal Differential Analysis
}
\description{
  Principal differential analysis (PDA) estimates a system of \eqn{n}
  linear differential equations that define functions that fit the data
  and their derivatives.  There is an equation in the system for each
  variable.

  Each equation has on its right side the highest order derivative that
  is used, and the order of this derivative, \eqn{m_j, j=1,...,n} can
  vary over equations.

  On the left side of equation is a linear combination of all
  the  variables and all the derivatives of these variables up to order
  one less than the order \eqn{m_j} of the highest derivative.

  In addition, the right side may contain linear combinations of forcing
  functions as well, with the number of forcing functions varying over
  equations.

  The linear combinations are defined by weighting functions multiplying
  each variable, derivative, and forcing function in the equation.
  These weighting functions may be constant or vary over time.  They are
  each represented by a functional parameter object, specifying a basis
  for an expansion of a coefficient, a linear differential operator for
  smoothing purposes, a smoothing parameter value, and a logical
  variable indicating whether the function is to be estimated, or kept
  fixed.
}
\usage{
pda.fd(xfdlist, bwtlist=NULL,
       awtlist=NULL, ufdlist=NULL, nfine=501)
}
\arguments{
  \item{xfdlist}{
    a list whose members are functional data objects representing each
    variable in the system of differential equations.  Each of these
    objects contain one or more curves to be represented by the
    corresponding differential equation.  The length of the list is
    equal to the number of differential equations. The number \eqn{N} of
    replications must be the same for each member functional data
    object.
  }
  \item{bwtlist}{
    this argument contains the weight coefficients that multiply, in the
    right side of each equation, all the variables in the system, and
    all their derivatives, where the derivatives are used up to one less
    than the order of the variable.   This argument has, in general, a
    three-level structure, defined by a three-level hierarchy of list
    objects.

    At the top level, the argument is a single list of length equal to
    the number of variables. Each component of this list is itself a
    list

    At the second level, each component of the top level list is itself
    a list, also of length equal to the number of variables.

    At the third and bottom level, each component of a second level list
    is a list of length equal to the number of orders of derivatives
    appearing on the right side of the equation, including  the variable
    itself, a derivative of order 0.  If m indicates the order of the
    equation, that is the order of the derivative on the left side, then
    this list is length m.

    The components in the third level lists are functional parameter
    objects defining estimates for weight functions.  For a first order
    equation, for example, \eqn{m = 1} and the single component in each
    list contains a weight function for the variable.  Since each
    equation has a term involving each variable in the system, a system
    of first order equations will have \eqn{n^2} at the third level of
    this structure.

    There MUST be a component for each weight function, even if the
    corresponding term does not appear in the equation.  In the case of
    a missing term, the corresponding component can be \code{NULL}, and
    it will be treated as a coefficient fixed at 0.

    However, in the case of a single differential equation,
    \code{bwtlist} can be given a simpler structure, since in this case
    only \eqn{m} coefficients are required.  Therefore, for a single
    equation, \code{bwtlist} can be a list of length \eqn{m} with each
    component containing a functional parameter object for the
    corresponding derivative.
  }

  \item{awtlist}{
    a two-level list containing weight functions for forcing functions.

    In addition to terms in each of the equations involving terms
    corresponding to each derivative of each variable in the system,
    each equation can also have a contribution from one or more
    exogenous variables, often called \emph{forcing functions.}

    This argument defines the weights multiplying these forcing
    functions, and is a list of length \eqn{n}, the number of
    variables.  Each component of this is is in turn a list, each
    component of which contains a functional parameter object defining a
    weight function for a forcing function.  If there are no forcing
    functions for an equation, this list can be \code{NULL}.  If none of
    the equations involve forcing functions, \code{awtlist} can be
    \code{NULL}, which is its default value if it is not in the argument
    list.
  }

  \item{ufdlist}{
    a two-level list containing forcing functions.  This list structure
    is identical to that for \code{awtlist}, the only difference being
    that the components in the second level contain functional data
    objects for the forcing functions, rather than functional parameter
    objects.
  }

  \item{nfine}{
    a number of values for a fine mesh.  The estimation of the
    differential equation involves discrete numerical quadrature
    estimates of integrals, and these require that functions be
    evaluated at a fine mesh of values of the argument.  This argument
    defines the number to use.  The default value of 501 is reset to
    five times the largest number of basis functions used to represent
    any variable in the system, if this number is larger.
  }
}

\value{
  an object of class \code{pda.fd}, being a list with the following
  components:

  \item{bwtlist}{
    a list array of the same dimensions as the corresponding argument,
    containing the estimated or fixed weight functions defining the
    system of linear differential equations.
  }
  \item{resfdlist}{
    a list of length equal to the number of variables or equations.
    Each members is a functional data object giving the residual
    functions or forcing functions defined as the left side of the
    equation (the derivative of order m of a variable) minus the linear
    fit on the right side.  The number of replicates for each residual
    functional data object is the same as that for the variables.
  }
  \item{awtlist}{
    a list of the same dimensions as the corresponding argument.  Each
    member is an estimated or fixed weighting function for a forcing
    function.
  }
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{pca.fd}},
  \code{\link{cca.fd}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
#See analyses of daily weather data for examples.
##
##  set up objects for examples
##

#  set up basis objects
#  constant basis object  for estimating weight functions
cbasis = create.constant.basis(c(0,1))
#  monomial basis: {1,t}  for estimating weight functions
mbasis = create.monomial.basis(c(0,1),2)
#  quartic spline basis with 54 basis functions for
#    defining functions to be analyzed
xbasis = create.bspline.basis(c(0,1),24,5)
#  set up functional parameter objects for weight bases
cfd0   = fd(0,cbasis)
cfdPar = fdPar(cfd0)
mfd0   = fd(matrix(0,2,1),mbasis)
mfdPar = fdPar(mfd0)

#  fine mesh for plotting functions
#  sampling points over [0,1]
tvec = seq(0,1,len=101)

##
##  Example 1:  a single first order constant coefficient unforced equation
##     Dx = -4*x  for  x(t) = exp(-4t)

beta    = 4
xvec    = exp(-beta*tvec)
xfd     = smooth.basis(tvec, xvec, xbasis)$fd
xfdlist = list(xfd)
bwtlist = list(cfdPar)
#  perform the principal differential analysis
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
##
##  Example 2:  a single first order varying coefficient unforced equation
##     Dx(t) = -t*x(t) or x(t) = exp(-t^2/2)
bvec    = tvec
xvec    = exp(-tvec^2/2)
xfd     = smooth.basis(tvec, xvec, xbasis)$fd
xfdlist = list(xfd)
bwtlist = list(mfdPar)
#  perform the principal differential analysis
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
##
##  Example 3:  a single second order constant coefficient unforced equation
##     Dx(t) = -(2*pi)^2*x(t) or x(t) = sin(2*pi*t)
##
xvec    = sin(2*pi*tvec)
xfd     = smooth.basis(tvec, xvec, xbasis)$fd
xfdlist = list(xfd)
bwtlist = list(cfdPar,cfdPar)
#  perform the principal differential analysis
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
##
##  Example 4:  two first order constant coefficient unforced equations
##     Dx1(t) = x2(t) and Dx2(t) = -x1(t)
##   equivalent to  x1(t) = sin(2*pi*t)
##
xvec1     = sin(2*pi*tvec)
xvec2     = 2*pi*cos(2*pi*tvec)
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
#  perform the principal differential analysis
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
##
##  Example 5:  a single first order constant coefficient equation
##     Dx = -4*x  for  x(t) = exp(-4t) forced by u(t) = 2
##
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
#  perform the principal differential analysis
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
##
##  Example 6:  two first order constant coefficient equations
##     Dx = -4*x    for  x(t) = exp(-4t)     forced by u(t) =  2
##     Dx = -4*t*x  for  x(t) = exp(-4t^2/2) forced by u(t) = -1
##
beta    = 4
xvec10  = exp(-beta*tvec)
alpha1  = 2
alpha2  = -1
xvec1   = xvec0 + alpha1*(1-xvec10)/beta
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

#  perform the principal differential analysis
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
awtfd1     = awtlistout[[1]][[1]]
awtfd2     = awtlistout[[2]][[1]]
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
par(oldpar)
}
\keyword{smooth}
