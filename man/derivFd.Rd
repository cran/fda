\name{derivFd}
\alias{derivFd}
\title{Compute a Derivative Functional Data Object }
\description{
A derivative of a functional data object, or the result of applying a linear differential operator to a functional data object, is then converted to a functional data object. This is intended for situations where a derivative is to be manipulated as a functional data object rather than simply evaluated. 
}
\usage{
derivFd(fd, Lfd=1)
}
\arguments{
  \item{fd }{A functional data object. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
}
\details{
The derivative of fd is first evaluated using function eval.fd, and then these values are converted to a functional data object using the data2fd function and the basis used to represent the original function. This assumes that the derivative does not vary any more complexly than the original function, which sometimes is doubtful. If the derivative is expected to be highly variable, it might be better to do this manually using a richer basis for the derivative rather than this function. 
}
\value{
A functional data object that represents the derivative of the argument fd. 
}
\references{}
\note{}
\author{}
\seealso{
eval.fd 
}
\examples{}
\keyword{}
