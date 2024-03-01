\name{smooth.surp}
\alias{smooth.surp}
\title{
  Fit data with surprisal smoothing.
}
\description{
  Surprisal is -log(probability) where the logarithm is to the base being the dimension 
  \code{M}of the multinomial observation vector.  The surprisal curves for each question are estimated
  by fitting the surprisal values of binned data using curves whose values are within the
  M-1 dimensional surprisal subspace that is within the space of non-negative M-dimensional vectors.
}
\usage{
smooth.surp(argvals, y, Bmat0, WfdPar, wtvec=NULL, conv=1e-4,
            iterlim=50, dbglev=0)
}
\arguments{
  \item{argvals}{
    Argument value array of length N, where N is the number of observed
    curve values for each curve.  It is assumed that that these argument
    values are common to all observed curves.  If this is not the case,
    you will need to run this function inside one or more loops,
    smoothing each curve separately.
  }
  \item{y}{
    A \code{nbin} by \code{M_i} matrix of surprisal values to be fit.  
  }
  \item{Bmat0}{A \code{Wnbasis} by \code{M_i - 1} matrix containing starting
  values for the iterative optimization of the least squares fit of the 
  surprisal curves to the surprisal data.
  }
  \item{WfdPar}{
    A functional parameter or fdPar object.  This object contains the
    specifications for the functional data object to be estimated by
    smoothing the data.  See comment lines in function fdPar for
    details.  The functional data object WFD in WFDPAROBJ is used to
    initialize the optimization process.  Its coefficient array contains
    the starting values for the iterative minimization of mean squared
    error.
  }
  \item{wtvec}{
    A vector of weights to be used in the smoothing.
  }
  \item{conv}{
    A convergence criterion.
  }
  \item{iterlim}{
    the maximum number of iterations allowed in the minimization of
    error sum of squares.
  }
  \item{dbglev}{
    Either 0, 1, or 2.  This controls the amount information printed out
    on each iteration, with 0 implying no output, 1 intermediate output
    level, and 2 full output.  If either level 1 or 2 is specified, it
    can be helpful to turn off the output buffering feature of S-PLUS.
  }
}
\value{
  A named list \code{surpFd} with these members:

  \item{Wfdobj}{
    a functional data object defining function $W(x)$ that that
    optimizes the fit to the data of the positive function that it
    defines.
  }
  \item{Flist}{
    a named list containing three results for the final converged
    solution: (1) \bold{f}: the optimal function value being minimized,
    (2) \bold{grad}: the gradient vector at the optimal solution, and
    (3) \bold{norm}: the norm of the gradient vector at the optimal
    solution.
  }
  \item{argvals}{
    the corresponding input arguments
  }
  \item{y}{
    the corresponding input arguments
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
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{eval.surp}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
#  evaluation points
x = seq(-2,2,len=11)
#  evaluate a standard normal distribution function
p = pnorm(x)
#  combine with 1-p
mnormp = cbind(p,1-p)
#  convert to surprisal values
mnorms = -log2(mnormp)
#  plot the surprisal values
matplot(x, mnorms, type="l", lty=c(1,1), col=c(1,1), 
        ylab="Surprisal (2-bits)")
# add some log-normal error
mnormdata = exp(log(mnorms) + rnorm(11)*0.1)
#  set up a b-spline basis object
nbasis = 7
sbasis = create.bspline.basis(c(-2,2),nbasis)
#  define an initial coefficient matrix
cmat = matrix(0,7,1)
#  set up a fdPar object for suprisal smoothing
sfd = fd(cmat, sbasis)
sfdPar = fdPar(sfd, Lfd=2, lambda=0)
#  smooth the noisy data
result = smooth.surp(x, mnormdata, cmat, sfdPar)
#  plot the data and the fits of the two surprisal curves
xfine = seq(-2,2,len=51)
sfine = eval.surp(xfine, result$Wfd)
matplot(xfine, sfine, type="l", lty=c(1,1), col=c(1,1))
points(x, mnormdata[,1])
points(x, mnormdata[,2])
#  convert the surprisal fit values to probabilities
pfine = 2^(-sfine)
#  check that they sum to one
apply(pfine,1,sum)
par(oldpar)
}

