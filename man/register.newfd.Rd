\name{register.newfd}
\alias{register.newfd}
\title{
  Register Functional Data Objects with Pre-Computed Warping Functions
}
\description{ This function registers a new functional data object to pre-computed
warping functions.
}
\usage{
register.newfd(yfd, Wfd,type=c('direct','monotone','periodic'))
}
\arguments{
  \item{yfd}{
    a multivariate functional data object defining the functions to be
    registered with \code{Wfd}.
  }
  \item{Wfd}{
    a functional data object defining the registration functions to be used
    to register \code{yfd}. This can be the result of either \code{landmarkreg}
    or \code{register.fd}.
  }
  \item{type}{
    indicates the type of registration function.
    \describe{
      \item{direct}{
	assumes \code{Wfd} is a direct definition of the registration
	functions. This is produced by \code{landmarkreg}.
      }
      \item{monotone}{
	assumes that \code{Wfd} defines a monotone functional data
	objected, up to shifting and scaling to make endpoints
	agree. This is produced by \code{register.fd}.
      }
      \item{periodic}{
	does shift registration for periodic functions. This is output
	from \code{register.fd} if \code{periodic=TRUE}.
      }
    }
  }
}
\value{
  a functional data object defining the registered curves.
}
\details{
  Only shift registration is considered for the periodic case.
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
  \code{\link{landmarkreg}},
  \code{\link{register.fd}}
}
