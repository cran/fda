\name{argvalsy.swap}
\alias{argvalsy.swap}
\title{
  Swap argvals with y if the latter is simpler.  
}
\description{
  Preprocess \code{argvals}, \code{y}, and \code{basisobj}.  If only one
  of \code{argvals} and \code{y} is provided, use  it as \code{y} and
  take \code{argvals} as a vector spanning basisobj[['rangreval']].  If
  both are provided, the simpler becomes \code{argvals}.  If both have
  the same dimensions but only one lies in basisobj[['rangreval']], that
  becomes \code{argvals}. 
}
\usage{
argvalsy.swap(argvals=NULL, y=NULL, basisobj=NULL) 
}
\arguments{
  \item{argvals}{
    a vector or array of argument values.
  }
  \item{y}{
    an array containing sampled values of curves.  
  }
  \item{basisobj}{
    One of the following:

    \itemize{
      \item{basisfd}{
	a functional basis object (class \code{basisfd}. 
      }
      \item{fd}{
	a functional data object (class \code{fd}), from which its
	\code{basis} component is extracted.  
      }
      \item{fdPar}{
	a functional parameter object (class \code{fdPar}), from which
	its \code{basis} component is extracted.  
      }
      \item{integer}{
	an integer giving the order of a B-spline basis,
	create.bspline.basis(argvals, norder=basisobj) 
      }
      \item{numeric vector}{
	specifying the knots for a B-spline basis,
	create.bspline.basis(basisobj) 
      } 	
      \item{NULL}{
	Defaults to create.bspline.basis(argvals).
      }
    }
  }
}
\details{
  1.  If \code{y} is NULL, replace by \code{argvals}.

  2.  If \code{argvals} is NULL, replace by
  seq(basisobj[['rangeval']][1], basisobj[['rangeval']][2], dim(y)[1])
  with a warning.  

  3.  If the dimensions of \code{argvals} and \code{y} match and only
  one is contained in basisobj[['rangeval']], use that as \code{argvals}
  and the other as \code{y}.

  4.  if \code{y} has fewer dimensions than \code{argvals}, swap them.  
}
\value{
  a list with components \code{argvals}, \code{y}, and \code{basisobj}.
}
\seealso{
  \code{\link{Data2fd}} 
  \code{\link{smooth.basis}}, 
  \code{\link{smooth.basisPar}}
}
\examples{
##
## one argument:  y
##
argvalsy.swap(1:5)
# warning ... 

##
## (argvals, y), same dimensions:  retain order 
##
argy1 <- argvalsy.swap(seq(0, 1, .2), 1:6)
argy1a <- argvalsy.swap(1:6, seq(0, 1, .2))

\dontshow{stopifnot(}
all.equal(argy1[[1]], argy1a[[2]]) &&
all.equal(argy1[[2]], argy1a[[1]])
# TRUE;  basisobj different 
\dontshow{)}

# lengths do not match 
\dontrun{
argvalsy.swap(1:4, 1:5)
} 

##
## two numeric arguments, different dimensions:  put simplest first 
##
argy2 <- argvalsy.swap(seq(0, 1, .2), matrix(1:12, 6))

\dontshow{stopifnot(}
all.equal(argy2,
argvalsy.swap(matrix(1:12, 6), seq(0, 1, .2)) )
# TRUE with a warning ... 
\dontshow{)}

\dontrun{
argvalsy.swap(seq(0, 1, .2), matrix(1:12, 2))
# ERROR:  first dimension does not match 
}

##
## one numeric, one basisobj
##
argy3 <- argvalsy.swap(1:6, b=4)
# warning:  argvals assumed seq(0, 1, .2) 

argy3. <- argvalsy.swap(1:6, b=create.bspline.basis(breaks=0:1))
# warning:  argvals assumed seq(0, 1, .2) 

argy3.6 <- argvalsy.swap(seq(0, 1, .2), b=create.bspline.basis(breaks=1:3))
# warning:  argvals assumed seq(1, 3 length=6)

##
## two numeric, one basisobj:  first matches basisobj
##
#  OK 
argy3a <- argvalsy.swap(1:6, seq(0, 1, .2),
              create.bspline.basis(breaks=c(1, 4, 8))) 

#  Swap (argvals, y) 
\dontshow{stopifnot(}
all.equal(argy3a,
argvalsy.swap(seq(0, 1, .2), 1:6, 
              create.bspline.basis(breaks=c(1, 4, 8))) )
# TRUE with a warning 
\dontshow{)}


\dontrun{
# neither match basisobj:  error  
argvalsy.swap(seq(0, 1, .2), 1:6, 
              create.bspline.basis(breaks=1:3) ) 
}

}
% docclass is function
\keyword{smooth}
