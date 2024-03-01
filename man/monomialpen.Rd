\name{monomialpen}
\alias{monomialpen}
\title{
  Evaluate Monomial Roughness Penalty Matrix
}
\description{
The roughness penalty matrix is the set of
inner products of all pairs of a derivative of integer powers of the
argument.
}
\usage{
monomialpen(basisobj, Lfdobj=int2Lfd(2),
            rng=basisobj$rangeval)
}
\arguments{
\item{basisobj}{
a monomial basis object.
}
\item{Lfdobj}{
either a nonnegative integer specifying an order of derivative
or a linear differential operator object.
}
\item{rng}{
the inner product may be computed over a range that is contained
within the range defined in the basis object.  This is a vector
or length two defining the range.
}
}
\value{
a symmetric matrix of order equal to the number of
monomial basis functions.
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
  \code{\link{exponpen}},
  \code{\link{fourierpen}},
  \code{\link{bsplinepen}},
  \code{\link{polygpen}}
}
\examples{
##
## set up a monomial basis for the first five powers
##
nbasis   <- 5
basisobj <- create.monomial.basis(c(-1,1),nbasis)
#  evaluate the rougness penalty matrix for the
#  second derivative.
penmat <- monomialpen(basisobj, 2)

##
## with rng of class Date and POSIXct
##
# Date
invasion1 <- as.Date('1775-09-04')
invasion2 <- as.Date('1812-07-12')
earlyUS.Canada <- c(invasion1, invasion2)
BspInvade1 <- create.monomial.basis(earlyUS.Canada)
invadmat <- monomialpen(BspInvade1)

# POSIXct
AmRev.ct <- as.POSIXct1970(c('1776-07-04', '1789-04-30'))
BspRev1.ct <- create.monomial.basis(AmRev.ct)
revmat <- monomialpen(BspRev1.ct)

}
% docclass is function
\keyword{smooth}
