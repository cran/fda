\name{symsolve}
\alias{symsolve}
\title{
  solve(A, B) where A is symmetric
}
\description{
  Solve A X = B for X where A is symmetric
}
\usage{
symsolve(Asym, Bmat)
}
\arguments{
  \item{Asym}{ a symmetric matrix }
  \item{Bmat}{ a square matrix of dimensions matching Asym }
}
\value{
  A square matrix of the same dimenions as Asym and Bmat
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
  \code{\link{solve}}
}
\examples{
A <- matrix(c(2,1,1,2), 2)
Ainv <- symsolve(A, diag(2))
}
% docclass is function
\keyword{array}
