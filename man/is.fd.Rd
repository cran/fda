\name{is.fd}
\alias{is.fd}
\title{
  Confirm Object has Class "fd"
}
\description{
Check that an argument is a functional data object.
}
\usage{
is.fd(fdobj)
}
\arguments{
\item{fdobj}{
an object to be checked.
}
}
\value{
a logical value:
\code{TRUE} if the class is correct, \code{FALSE} otherwise.
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
  \code{\link{is.basis}}, 
  \code{\link{is.fdPar}}, 
  \code{\link{is.Lfd}}
}
% docclass is function
\keyword{smooth}
