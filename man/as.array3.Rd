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
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
  \emph{Functional data analysis with R and Matlab}, Springer, New
  York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
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
  as.array3(matrix(1:6, 2, dimnames=list(letters[1:2],
      LETTERS[3:5]))) 

  ##
  ## array -> array 
  ##
  
  as.array3(array(1:6, 1:3)) 

  ##
  ## 4-d array 
  ##
  # These lines throw an error because the dimensionality woud be 4
  # and as.array3 only allows dimensions 3 or less.
  # if(!CRAN()) {
  #   as.array3(array(1:24, 1:4)) 
  # }
}
\keyword{utilities}

