\name{create.fourier.basis}
\alias{create.fourier.basis}
\title{Create a set of Fourier basis functions for constructing a functional data object. }
\description{
Functional data objects are constructed by specifying a set of basis functions and a set of coefficients defining a linear combination of these basis functions. The Fourier basis is a system that is usually used for periodic functions. It has the advantages of very fast computation and great flexibility. If the data are considered to be nonperiod, the Fourier basis is usually preferred. 
}
\usage{
create.fourier.basis(rangeval=c(0, 1), nbasis=3, period=width)
}
\arguments{
  \item{rangeval}{A vector of length 2 containing the initial and final values of the interval over which the functional data object can be evaluated. Default value c(0,1) }
  \item{nbasis }{An integer variable specifying the number of basis functions. Default value 3. }
  \item{period }{The width of any interval over which the Fourier functions repeat themselves, or are periodic. The default is the width of the interval defined in rangeval. }
}
\details{
The first basis function is the unit function with the value one everywhere. The next two are the sine/cosine pair with period defined in the argument period. The fourth and fifth are the sin/cosine series with period one half of period. And so forth. The number of basis functions is usually odd. 
}
\value{
Returns a list with the basis.fd, class attribute with entries as above having names type, rangeval, nbasis, and params, respectively. The type member is fourier. 
}
\references{}
\note{}
\author{}
\seealso{
basis.fd. create.basis, create.bspline.basis 
}
\examples{}
\keyword{}
