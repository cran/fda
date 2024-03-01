\name{zerofind}
\alias{zerofind}
\title{
  Does the range of the input contain 0?  
}
\description{
  Returns TRUE if range of the argument includes 0 and FALSE if not.  
}
\usage{
  zerofind(fmat)
}
\arguments{
  \item{fmat}{
    An object from which 'range' returns two numbers.  
  }
}
\value{
  A logical value TRUE or FALSE.  
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
  \code{\link{range}}
}
\examples{
zerofind(1:5)
# FALSE
zerofind(0:3)
# TRUE 
}
% docclass is function
\keyword{logic}
