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
\seealso{
  \code{\link{solve}}
}
\examples{
A <- matrix(c(2,1,1,2), 2)
Ainv <- symsolve(A, diag(2))
}
% docclass is function
\keyword{array}
