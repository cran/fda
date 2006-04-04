\name{create.constant.basis}
\alias{create.constant.basis}
\title{
  Create a Constant Basis
}
\description{
Create a constant basis object, defining a single basis function
whose value is everywhere 1.0.
}
\usage{
create.constant.basis(rangeval=c(0, 1))
}
\arguments{
\item{rangeval}{
a vector of length 2 containing the initial and final
values of argument t defining the interval over which the functional
data object can be evaluated.
However, this is seldom used
since the value of the basis function does not depend on the range or
any argument values.
}
}
\value{
a basis object of the type \code{const}.
}
\seealso{
\code{\link{basisfd}}, 
\code{\link{create.bspline.basis}}, 
\code{\link{create.exponential.basis}}, 
\code{\link{create.fourier.basis}}, 
\code{\link{create.monomial.basis}}, 
\code{\link{create.polygonal.basis}}, 
\code{\link{create.polynomial.basis}}, 
\code{\link{create.power.basis}}
}
\examples{

basisobj <- create.constant.basis(c(-1,1))

}
% docclass is function
\keyword{smooth}
