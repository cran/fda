\name{smooth.fd}
\alias{smooth.fd}
\title{Smooth a Functional Data Object Using a Roughness Penalty }
\description{
This function is intended to apply a roughness penalized smooth to data already set up as a functional data object. For example, data may have been converted to a functional data object using function data2fd using a fairly large set of basis functions, and subsequently it was desired to smooth the functional data object that resulted. 
}
\usage{
smooth.fd(fd, lambda=0, Lfd=NULL, rebase=TRUE)
}
\arguments{
  \item{fd }{The functional data object to be smoothed. }
  \item{lambda }{A nonnegative value controlling the amount of roughness in the data. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
  \item{rebase }{A logical variable that is only relevant if the basis is a polygonal basis of type "polyg". If this case, if rebase is TRUE , then the basis is changed to a cubic bspline basis before smoothing. }
}
\details{
See function smooth.basis for details. 
}
\value{
A functional data object. 
}
\references{
A discussion of roughness penalties can be found in Chapter 4 of Ramsay, J. O. and Silverman, B.W. (1997) Functional Data Analysis. More information can be found in recent texts on nonparametric regression. 
}
\note{}
\author{}
\seealso{
data2fd, plotfit.fd, smooth.basis, project.basis 
}
\examples{}
\keyword{}
