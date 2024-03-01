\name{sparse.list}
\alias{sparse.list}
\title{Creates a list of sparse data from a matrix}
\description{
  Creates a list with sparse data from a matrix that has NAs. The length of 
  the list will be equal to the number of columns in the data matrix (replications)
}
\usage{
  sparse.list(data,time)
}
\arguments{
  \item{data}{
    If the set is supplied as a matrix object, the rows must correspond to argument
    values and columns to replications, and it will be assumed that there is only 
    one variable per observation.  If y is a three-dimensional array, the first
    dimension corresponds to argument values, the second to replications, and the
    third to variables within replications.
  }
  \item{time}{Time points where the observations where taken.}
}
\value{
  a list with \code{length dim(data)[2]}. Each element of the list is a matrix with \code{ncol > 1}. The first column of each element corresponds to the point index per observation.
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
  \code{\link{cca.fd}},
  \code{\link{pda.fd}}
}
