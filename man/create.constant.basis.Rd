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
basisfd, create.bspline.basis, create.exponential.basis, create.fourier.basis,
create.monomial.basis, create.polygonal.basis, create.polynomial.basis,
create.power.basis
}
\examples{

basisobj <- create.constant.basis(c(-1,1))

}
% docclass is function
\keyword{smooth}
