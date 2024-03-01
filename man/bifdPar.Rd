\name{bifdPar}
\alias{bifdPar}
\title{
  Define a Bivariate Functional Parameter Object
}
\description{
  Functional parameter objects are used as arguments to functions that
  estimate functional parameters, such as smoothing functions like
  \code{smooth.basis}.  A bivariate functional parameter object supplies
  the analogous information required for smoothing bivariate data using
  a bivariate functional data object $x(s,t)$.  The arguments are the same as
  those for \code{fdPar} objects, except that two linear differential
  operator objects and two smoothing parameters must be applied,
  each pair corresponding to one of the arguments $s$ and $t$ of the
  bivariate functional data object.
}
\usage{
bifdPar(bifdobj, Lfdobjs=int2Lfd(2), Lfdobjt=int2Lfd(2), lambdas=0, lambdat=0,
      estimate=TRUE)
}
\arguments{
  \item{bifdobj}{
    a bivariate functional data object.
  }
  \item{Lfdobjs}{
    either a nonnegative integer or a linear differential operator
    object for the first argument $s$.

    If \code{NULL}, Lfdobjs depends on bifdobj[['sbasis']][['type']]:

    \describe{
      \item{bspline}{
	Lfdobjs <- int2Lfd(max(0, norder-2)), where norder =
	norder(bifdobj[['sbasis']]).
      }
      \item{fourier}{
	Lfdobjs = a harmonic acceleration operator:

	\code{Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rngs))^2,0), rngs)}

	where rngs = bifdobj[['sbasis']][['rangeval']].
      }
      \item{anything else}{Lfdobj <- int2Lfd(0)}
    }
  }
  \item{Lfdobjt}{
    either a nonnegative integer or a linear differential operator
    object for the first argument $t$.

    If \code{NULL}, Lfdobjt depends on bifdobj[['tbasis']][['type']]:

    \describe{
      \item{bspline}{
	Lfdobj <- int2Lfd(max(0, norder-2)), where norder =
	norder(bifdobj[['tbasis']]).
      }
      \item{fourier}{
	Lfdobj = a harmonic acceleration operator:

	\code{Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rngt))^2,0), rngt)}

	where rngt = bifdobj[['tbasis']][['rangeval']].
      }
      \item{anything else}{Lfdobj <- int2Lfd(0)}
    }
  }
  \item{lambdas}{
    a nonnegative real number specifying the amount of smoothing
    to be applied to the estimated functional parameter $x(s,t)$
    as a function of $s$..
  }
  \item{lambdat}{
    a nonnegative real number specifying the amount of smoothing
    to be applied to the estimated functional parameter $x(s,t)$
    as a function of $t$..
  }
  \item{estimate}{not currently used.}
}
\value{
  a bivariate functional parameter object (i.e., an object of class 
  \code{bifdPar}), which is a list with the following components:

  \item{bifd}{
    a functional data object (i.e., with class \code{bifd})
  }
  \item{Lfdobjs}{
    a linear differential operator object (i.e., with class
    \code{Lfdobjs})
  }
  \item{Lfdobjt}{
    a linear differential operator object (i.e., with class
    \code{Lfdobjt})
  }
  \item{lambdas}{
    a nonnegative real number
  }
  \item{lambdat}{
    a nonnegative real number
  }
  \item{estimate}{not currently used}
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
\seealso{
\code{\link{linmod}}
}
\source{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009)
    \emph{Functional Data Analysis in R and Matlab}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York
}
\examples{
#See the prediction of precipitation using temperature as
#the independent variable in the analysis of the daily weather
#data, and the analysis of the Swedish mortality data.
}


% docclass is function
\keyword{bivariate smooth}
