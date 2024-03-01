\name{powerpen}
\alias{powerpen}
\title{
  Power Penalty Matrix
}
\description{
Computes the matrix defining the roughness penalty for functions
expressed in terms of a power basis.
}
\usage{
powerpen(basisobj, Lfdobj=int2Lfd(2))
}
\arguments{
\item{basisobj}{
a power basis object.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
}
}
\value{
a symmetric matrix of order equal to the number of basis functions
defined by the power basis object.  Each element is the inner product
of two power basis functions after applying the derivative or linear
differential operator defined by \code{Lfdobj}.
}
\details{
A roughness penalty for a function $ x(t) $ is defined by
integrating the square of either the derivative of  $ x(t) $ or,
more generally, the result of applying a linear differential operator
$ L $ to it.  The most common roughness penalty is the integral of
the square of the second derivative, and
this is the default. To apply this roughness penalty, the matrix of
inner products produced by this function is necessary.
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
  \code{\link{create.power.basis}}, 
  \code{\link{powerbasis}}
}
\examples{

#  set up an power basis with 3 basis functions.
#  the powers are 0, 1, and 2.
basisobj <- create.power.basis(c(0,1),3,c(0,1,2))
#  compute the 3 by 3 matrix of inner products of second derivatives
#FIXME
#penmat <- powerpen(basisobj, 2)

}
% docclass is function
\keyword{smooth}
