\name{as.array3}
\alias{as.array3}
\alias{as.array3}
\title{
  Reshape a vector or array to have 3 dimensions.  
}
\description{
  Coerce a vector or array to have 3 dimensions, preserving dimnames if
  feasible.  Throw an error if length(dim(x)) > 3.  
}
\usage{
as.array3(x) 
}
\arguments{
  \item{x}{
    A vector or array.  
  }
}
\details{
  1.  dimx <- dim(x);  ndim <- length(dimx) 

  2.  if(ndim==3)return(x).

  3.  if(ndim>3)stop.

  4.  x2 <- as.matrix(x)

  5.  dim(x2) <- c(dim(x2), 1)

  6.  xnames <- dimnames(x)

  7.  if(is.list(xnames))dimnames(x2) <- list(xnames[[1]], xnames[[2]],
  NULL) 
}
\value{
  A 3-dimensional array with names matching \code{x} 
}
%\references{ }
\author{
  Spencer Graves
} 
\seealso{
  \code{\link{dim}},
  \code{\link{dimnames}}
  \code{\link{checkDims3}}
}
\examples{
##
## vector -> array 
##
as.array3(c(a=1, b=2)) 

##
## matrix -> array 
##
as.array3(matrix(1:6, 2))
as.array3(matrix(1:6, 2, dimnames=list(letters[1:2], LETTERS[3:5]))) 

##
## array -> array 
##
as.array3(array(1:6, 1:3)) 

##
## 4-d array 
##
\dontrun{
as.array3(array(1:24, 1:4)) 
Error in as.array3(array(1:24, 1:4)) : 
  length(dim(array(1:24, 1:4)) = 4 > 3
}
}
\keyword{utilities}

