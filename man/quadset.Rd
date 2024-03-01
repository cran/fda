\name{quadset}
\alias{quadset}
\title{
  Quadrature points and weights for Simpson's rule
}
\description{
  Set up quadrature points and weights for Simpson's rule.
}
\usage{
quadset(nquad=5, basisobj=NULL, breaks)
}
\arguments{
  \item{nquad}{
    an odd integer at least 5 giving the number of evenly spaced
    Simpson's rule quadrature points to use over each interval
    (breaks[i], breaks[i+1]).
  }
  \item{basisobj}{
    the basis object that will contain the quadrature points and weights
  }
  \item{breaks}{
    optional interval boundaries.  If this is provided, the first value
    must be the initial point of the interval over which the basis is
    defined, and the final value must be the end point.  If this is not
    supplied, and 'basisobj' is of type 'bspline', the knots are used as
    these values.
  }
}
\details{
  Set up quadrature points and weights for Simpson's rule and store
  information in the output 'basisobj'.  Simpson's rule is used to
  integrate a function between successive values in vector 'breaks'.
  That is, over each interval (breaks[i], breaks[i+1]).  Simpson's rule
  uses 'nquad' equally spaced quadrature points over this interval,
  starting with the the left boundary and ending with the right
  boundary.  The quadrature weights are the values
  delta*c(1,4,2,4,2,4,..., 2,4,1) where 'delta' is the difference
  between successive quadrature points, that is, delta =
  (breaks[i-1]-breaks[i])/(nquad-1).
}
\value{
  If is.null(basisobj), quadset returns a 'quadvals' matrix with columns
  quadpts and quadwts.  Otherwise, it returns basisobj with the
  two components set as follows:

  \item{quadvals}{
    cbind(quadpts=quadpts, quadwts=quadwts)
  }
  \item{value}{
    a list with two components containing eval.basis(quadpts, basisobj,
    ival-1) for ival=1, 2.
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
  \code{\link{create.bspline.basis}}
  \code{\link{eval.basis}}
}
\examples{
(qs7.1 <- quadset(nquad=7, breaks=c(0, .3, 1)))
# cbind(quadpts=c(seq(0, 0.3, length=7),
#              seq(0.3, 1, length=7)),
#    quadwts=c((0.3/18)*c(1, 4, 2, 4, 2, 4, 1),
#              (0.7/18)*c(1, 4, 2, 4, 2, 4, 1) ) )

# The simplest basis currently available with this function:
bspl2.2 <- create.bspline.basis(norder=2, breaks=c(0,.5, 1))

bspl2.2a <- quadset(basisobj=bspl2.2)
bspl2.2a$quadvals
# cbind(quadpts=c((0:4)/8, .5+(0:4)/8),
#       quadwts=rep(c(1,4,2,4,1)/24, 2) )
bspl2.2a$values
# a list of length 2
# [[1]] = matrix of dimension c(10, 3) with the 3 basis
#      functions evaluated at the 10 quadrature points:
# values[[1]][, 1] = c(1, .75, .5, .25, rep(0, 6))
# values[[1]][, 2] = c(0, .25, .5, .75, 1, .75, .5, .25, 0)
# values[[1]][, 3] = values[10:1, 1]
#
# values[[2]] = matrix of dimension c(10, 3) with the
#     first derivative of values[[1]], being either
#    -2, 0, or 2.
}
% docclass is function
\keyword{smooth}
