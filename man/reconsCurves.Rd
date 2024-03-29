\name{reconsCurves}
\alias{reconsCurves}
\title{
  Reconstruct data curves using functional principal components 
}
\description{
  Reconstructs data curves as objects of class \code{fd} using functional principal   components
}
\usage{
  reconsCurves(data, PC)
}
\arguments{
  \item{data}{
    a set of values of curves at discrete sampling points or
    argument values. If the set is supplied as a matrix object, the rows must
    correspond to argument values and columns to replications, and it will be
    assumed that there is only one variable per observation.  If
    \code{data} is a three-dimensional array, the first dimension
    corresponds to argument values, the second to replications, and the
    third to variables within replications.
  }
  \item{PC}{
    an object of class \code{pca.fd}.
  }
}
\value{
  an object of class \code{fd}.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
