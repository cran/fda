\name{exponpen}
\alias{exponpen}
\title{
  Exponential Penalty Matrix
}
\description{
Computes the matrix defining the roughness penalty for functions
expressed in terms of an exponential basis.
}
\usage{
exponpen(basisobj, Lfdobj=int2Lfd(2))
}
\arguments{
\item{basisobj}{
an exponential basis object.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
}
}
\value{
a symmetric matrix of order equal to the number of basis functions
defined by the exponential basis object.  Each element is the inner product
of two exponential basis functions after applying the derivative or linear
differential operator defined by Lfdobj.
}
\details{
A roughness penalty for a function $x(t)$ is defined by
integrating the square of either the derivative of  $ x(t) $ or,
more generally, the result of applying a linear differential operator
$L$ to it.  The most common roughness penalty is the integral of
the square of the second derivative, and
this is the default. To apply this roughness penalty, the matrix of
inner products of the basis functions (possibly after applying the
linear differential operator to them) defining this function
is necessary. This function just calls the roughness penalty evaluation
function specific to the basis involved.
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
  \code{\link{expon}}, 
  \code{\link{eval.penalty}}, 
  \code{\link{getbasispenalty}}
}
\examples{

#  set up an exponential basis with 3 basis functions
ratevec  <- c(0, -1, -5)
basisobj <- create.exponential.basis(c(0,1),3,ratevec)
#  compute the 3 by 3 matrix of inner products of
#  second derivatives
penmat <- exponpen(basisobj)

}
% docclass is function
\keyword{smooth}
