\name{sparse.mat}
\alias{sparse.mat}
\title{Creates a matrix of sparse data with NAs out of a list}
\description{
  Creates a matrix or three dimensional array of sparse data with NAs from a list of
  sparse data. The number of columns of the matrix will be equal to the length of the
  list (replications)
}
\usage{
  sparse.mat(datalist)
}
\arguments{
  \item{datalist}{
    A list object. Each element must be a matrix object where the rows correspond to 
    argument values per individual. First column corresponds to time points and the 
    following columns to argument values per variable.
  }
}
\value{
  a matrix or array with ncol = length(datalist). First dimension corresponds to point
  observations, second dimension corresponds to replications and third dimension
  corresponds to variables.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
