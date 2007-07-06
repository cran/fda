\name{register.fd}
\alias{register.fd}
\title{
  Register Functional Data Objects Using a Continuous Criterion
}
\description{

criterion.  By aligned is meant that the shape of each curve is matched
as closely as possible to that of the target by means of a smooth
increasing transformation of the argument, or a warping function.
}
\usage{
register.fd(y0fd=NULL, yfd=NULL, WfdParobj=c(Lfdobj=2, lambda=1),
            conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)
}
\arguments{
  \item{y0fd}{
    a functional data object defining the target for registration.

    If \code{yfd} is NULL and y0fd is a multivariate data object, then
    y0fd is assigned to yfd and y0fd is replaced by its mean.

    Alternatively, if \code{yfd} is a multivariate functional data
    object and y0fd is missing, y0fd is replaced by the mean of
    \code{y0fd}.  

    Otherwise, y0fd must be a univariate functional data object taken as
    the target to which \code{yfd} is registered.  
  }
  \item{yfd}{
    a multivariate functional data object defining the functions to be
    registered to target \code{y0fd}.  If it is NULL and \code{y0fd} is
    a multivariate functional data object, yfd takes the value of 
    \code{y0fd}.  
  }
  \item{WfdParobj}{
    a functional parameter object for a single function.  This is used
    as the initial value in the estimation of a function $W(t)$ that
    defines the warping function $h(t)$ that registers a particular
    curve. The object also contains information on a roughness penalty 
    and smoothing parameter to control the roughness of $h(t)$.

    Alternatively, this can be a vector or a list with components named
    \code{Lfdobj} and \code{lambda}, which are passed as arguments to
    \code{fdPar} to create the functional parameter form of WfdParobj
    required by the rest of the register.fd algorithm.

    The default \code{Lfdobj} of 2 penalizes curvature, thereby
    preferring no warping of time, with \code{lambda} indicating the
    strength of that preference.  A common alternative is \code{Lfdobj}
    = 3, penalizing the rate of change of curvature.   
  }
  \item{conv}{
    a criterion for convergence of the iterations.
  }
  \item{iterlim}{
    a limit on the number of iterations.
  }
  \item{dbglev}{
    either 0, 1, or 2.  This controls the amount information printed out
    on each iteration, with 0 implying no output, 1 intermediate output
    level, and 2 full output.  (If this is run with output buffering
    such as used with S-Plus, it may be necessary to turn off the output
    buffering to actually get the progress reports before the completion
    of computations.)  
  }
  \item{periodic}{
    a logical variable:  if \code{TRUE}, the functions are considered to
    be periodic, in which case a constant can be added to all argument
    values after they are warped. 
  }
  \item{crit}{
    an integer that is either 1 or 2 that indicates the nature of the
    continuous registration criterion that is used.  If 1, the criterion is
    least squares, and if 2, the criterion is the minimum eigenvalue of a
    cross-product matrix.  In general, criterion 2 is to be preferred.
  }
}
\value{
  a named list of length 3 containing the following components:

  \item{regfd}{
    A functional data object containing the registered functions.
  }
  \item{Wfd}{
    A functional data object containing the functions $h W(t)$
    that define the warping functions $h(t)$.
  }
  \item{shift}{
    If the functions are periodic, this is a vector of time shifts.
  }
}
\details{
  The warping function that smoothly and monotonely transforms the
  argument is defined by \code{Wfd} is the same as that defines the
  monotone smoothing function in for function \code{smooth.monotone.}
  See the help file for that function for further details.
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York, ch. 6 & 7.  
}
\seealso{
  \code{\link{smooth.monotone}}, 
  \code{\link{smooth.morph}}
}
\examples{
#See the analyses of the growth data for examples.
##
## 1.  Simplest call
##
# Specify smoothing weight 
lambda.gr2.3 <- .03

# Specify what to smooth, namely the rate of change of curvature
Lfdobj.growth    <- 2 

# Establish a B-spline basis
nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
rng.growth <- range(growth$age)
# 1 18 
wbasis.growth <- create.bspline.basis(rangeval=rng.growth,
                   nbasis=nbasis.growth, norder=norder.growth,
                   breaks=growth$age)

# Smooth consistent with the analysis of these data
# in afda-ch06.R, and register to individual smooths:  
cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)
growfdPar2.3 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2.3)
# Create a functional data object for all the boys
hgtmfd.all <- with(growth, smooth.basis(age, hgtm, growfdPar2.3))

nBoys <- 2
# nBoys <- dim(growth[["hgtm"]])[2]
# register.fd takes time, so use only 2 curves as an illustration
# to minimize compute time in this example;  

#Alternative to subsetting later is to subset now:  
#hgtmfd.all<-with(growth,smooth.basis(age, hgtm[,1:nBoys],growfdPar2.3))

# Register the growth velocity rather than the
# growth curves directly 
smBv <- deriv(hgtmfd.all$fd, 1)

# This takes time, so limit the number of curves registered to nBoys

\dontrun{
smB.reg.0 <- register.fd(smBv[1:nBoys])

smB.reg.1 <- register.fd(smBv[1:nBoys],WfdParobj=c(Lfdobj=Lfdobj.growth, lambda=lambda.gr2.3))

##
## 2.  Call providing the target
##

smBv.mean <- deriv(mean(hgtmfd.all$fd[1:nBoys]), 1)
smB.reg.2a <- register.fd(smBv.mean, smBv[1:nBoys],
               WfdParobj=c(Lfdobj=Lfdobj.growth, lambda=lambda.gr2.3))

smBv.mean <- mean(smBv[1:nBoys]) 
smB.reg.2 <- register.fd(smBv.mean, smBv[1:nBoys],
               WfdParobj=c(Lfdobj=Lfdobj.growth, lambda=lambda.gr2.3))
all.equal(smB.reg.1, smB.reg.2) 

##
## 3.  Call using WfdParobj
##

# Create a dummy functional data object
# to hold the functional data objects for the
# time warping function
# ... start with a zero matrix (nbasis.growth, nBoys) 
smBc0 <- matrix(0, nbasis.growth, nBoys)
# ... convert to a functional data object 
smBwfd0 <- fd(smBc0, wbasis.growth)
# ... convert to a functional parameter object 
smB.wfdPar <- fdPar(smBwfd0, Lfdobj.growth, lambda.gr2.3)

smB.reg.3<- register.fd(smBv[1:nBoys], WfdParobj=smB.wfdPar)
all.equal(smB.reg.1, smB.reg.3)
}

}
% docclass is function
\keyword{smooth}
