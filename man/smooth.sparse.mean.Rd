\name{smooth.sparse.mean}
\alias{smooth.sparse.mean}

\title{
  Smooth the mean function of sparse data 
}
\description{
  Do a smoothing of the mean function for sparse data that is either given 
  as a list or as a matrix with NAs. The smooth is done by basis expansion 
  with the functional basis "type"; if !(lambda == 0) then the second 
  derivative is penalized (int2Lfd(2)).
}
\usage{
smooth.sparse.mean(data, time ,rng = c(0, 1), type = "", nbasis = NULL, 
                  knots = NULL, norder = NULL, lambda = NULL)
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
    column corresponds to time points and followins columns to argument values 
    per variable.
}
  \item{time}{
  Array with time points where data was taken. length(time) == ncol(data)
}
  \item{rng}{
  an array of length 2 containing the lower and upper boundaries for the 
  rangeval of argument values
}
  \item{type}{
  Type of basisfd for smoothing the mean estimate function.
  "bspline", "fourier", "exp", "const" or "mon"
}
  \item{nbasis}{
  An integer variable specifying the number of basis functions
}
  \item{knots}{
  a vector specifying the break points if type == "bspline"
}
  \item{norder}{
  an integer specifying the order of b-splines if type == "bspline"
}
  \item{lambda}{
  a nonnegative real number specifying the amount of smoothing to be applied to 
  the estimated functional parameter
  }
}
\value{
a functional data object containing a smooth of the mean.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}

