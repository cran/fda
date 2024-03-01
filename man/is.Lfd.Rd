\name{is.Lfd}
\alias{is.Lfd}
\title{
  Confirm Object has Class "Lfd"
}
\description{
Check that an argument is a linear differential operator object.
}
\usage{
is.Lfd(Lfdobj)
}
\arguments{
\item{Lfdobj}{
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
  \code{\link{is.fd}}, 
  \code{\link{is.fdPar}}
}
\keyword{smooth}
