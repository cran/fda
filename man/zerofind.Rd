\name{zerofind}
\alias{zerofind}
\title{
  Does the range of the input contain 0?  
}
\description{
  Returns TRUE if range of the argument includes 0 and FALSE if not.  
}
\usage{
  zerofind(fmat, na.rm=TRUE)
}
\arguments{
  \item{fmat}{
    An object from which 'range' returns two numbers.  
  }
  \item{na.rm}{
    logical:  TRUE to remove NAs before computing \code{\link{range}}
  }
}
\value{
  A logical value TRUE or FALSE.  
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
