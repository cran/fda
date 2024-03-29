\name{is.eqbasis}
\alias{is.eqbasis}
\title{
  Confirm that two objects of  class "Basisfd" are identical
}
\description{
Check all the slots of two basis objects to see that they are identical.
}
\usage{
is.eqbasis(basisobj1, basisobj2)
}
\arguments{
\item{basisobj1}{
The first basis object to be checked for being identical to the second.
}
\item{basisobj2}{
The second basis object to be checked for being identical to the first.
}
}
\value{
  a logical value:
  \code{TRUE} if the two basis objects are identical, \code{FALSE} otherwise.
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
  \code{\link{basisfd}},
}
% docclass is function
\keyword{basis}
