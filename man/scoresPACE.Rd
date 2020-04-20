\name{scoresPACE}
\alias{scoresPACE}
\title{Estimates of functional Principal Component scores through PACE}
\description{
  Function \code{scoresPACE} estimates the functional Principal Component 
  scores through Conditional Expectation (PACE)
}
\usage{
  scoresPACE(data, time, covestimate, PC)
}
\arguments{
  \item{data}{
    a matrix object or list -- If the set is supplied as a matrix object, 
    the rows must correspond to argument values and columns to replications, 
    and it will be assumed that there is only one variable per observation.  
    If y is a three-dimensional array, the first dimension corresponds to  
    argument values, the second to replications, and the third to variables 
    within replications. -- If it is a list, each element must be a matrix
    object, the rows correspond to argument values per individual. First 
    column corresponds to time points and following columns to argument values 
    per variable.
  }
  \item{time}{
    Array with time points where data was taken. \code{length(time) == dim(data)[1]}
  }
  \item{covestimate}{
    a list with the two named entries "cov.estimate" and "meanfd"
  }
  \item{PC}{
    an object of class "pca.fd"
  }
}
\value{
  a matrix of scores with dimension nrow = nharm and ncol = ncol(data)
}
