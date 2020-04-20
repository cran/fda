\name{covPACE}
\alias{covPACE}
\title{
  Estimate of the covariance surface 
}
\description{
  Function \code{covPACE} does a bivariate smoothing for estimating the 
  covariance surface for data that has not yet been smoothed
}
\usage{
  covPACE(data,rng , time, meanfd, basis, lambda, Lfdobj)
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
  \item{rng}{
    a vector of length 2 defining a restricted range where the data was observed
  }
  \item{time}{
    Array with time points where data was taken. \code{length(time) == dim(data)[1]}
  }
  \item{meanfd}{
    Fd object corresponding to the mean function of the data
  }
  \item{basis}{
    basisfd object for smoothing the covariate function
  }
  \item{lambda}{
    a nonnegative real number specifying the amount of smoothing to be applied to
    the estimated functional parameter
  }
  \item{Lfdobj}{
    a linear differential operator object for smoothing penalty of the estimated 
    functional parameter
  }
}
\value{
  a list with these two named entries:
    \item{cov.estimate}{
      an object of class "bifd" object or a list of "bifd" elements
    }
    \item{meanfd}{
      a functional data object giving the mean function
    }
}
