\name{checkLogicalInteger}
\alias{checkLogical}
\alias{checkNumeric}
\alias{checkLogicalInteger}
\title{
  Does an argument satisfy required conditions?    
}
\description{
  Check whether an argument is a logical vector of a certain length or a
  numeric vector in a certain range and issue an appropriate error or
  warning if not:

  \code{checkLogical} throws an error or returns FALSE with a warning
  unless \code{x} is a  logical vector of exactly the required
  \code{length}.
  
  \code{checkNumeric} throws an error or returns FALSE with a warning
  unless \code{x} is either NULL or a \code{numeric} vector of at most
  \code{length} with \code{x} in the desired range.  
  
  \code{checkLogicalInteger} returns a logical vector of exactly
  \code{length} unless \code{x} is neither NULL nor \code{logical} of
  the required \code{length} nor \code{numeric} with \code{x} in the
  desired range. 
}
\usage{
checkLogical(x, length., warnOnly=FALSE)
checkNumeric(x, lower, upper, length., integer=TRUE, unique=TRUE,
             inclusion=c(TRUE,TRUE), warnOnly=FALSE)
checkLogicalInteger(x, length., warnOnly=FALSE)
}
\arguments{
  \item{x}{ an object to be checked }
  \item{length.}{   
    The required length for \code{x} if \code{logical} and not NULL or
    the maximum length if \code{numeric}.  
  }
  \item{lower, upper}{
    lower and upper limits for \code{x}.  
  } 
  \item{integer}{
    logical:  If true, a \code{numeric} \code{x} must be
    \code{integer}.  
  }
  \item{unique}{
    logical:  TRUE if duplicates are NOT allowed in \code{x}.  
  }
  \item{inclusion}{
    logical vector of length 2, similar to
    \code{link[ifultools]{checkRange}}:  
    
    if(inclusion[1]) (lower <= x) else (lower < x)

    if(inclusion[2]) (x <= upper) else (x < upper)
  }   
  \item{warnOnly}{
    logical:  If TRUE, violations are reported as warnings, not as
    errors.  
  }
}
\details{
  1.  xName <- deparse(substitute(x)) to use in any required error or
  warning.  

  2.  if(is.null(x)) handle appropriately:  Return FALSE for
  \code{checkLogical}, TRUE for \code{checkNumeric} and rep(TRUE,
  length.) for \code{checkLogicalInteger}.  

  3.  Check class(x).

  4.  Check other conditions.  
}
\value{
  \code{checkLogical} returns a logical vector of the required
  \code{length.}, unless it issues an error message.

  \code{checkNumeric} returns a numeric vector of at most \code{length.}
  with all elements between \code{lower} and \code{upper}, and
  optionally \code{unique}, unless it issues an error message.

  \code{checkLogicalInteger} returns a logical vector of the required
  \code{length.}, unless it issues an error message.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\author{
  Spencer Graves
} 
\examples{
##
## checkLogical
##
checkLogical(NULL, length=3, warnOnly=TRUE)
checkLogical(c(FALSE, TRUE, TRUE), length=4, warnOnly=TRUE)
checkLogical(c(FALSE, TRUE, TRUE), length=3)

##
## checkNumeric
##
checkNumeric(NULL, lower=1, upper=3)
checkNumeric(1:3, 1, 3)
checkNumeric(1:3, 1, 3, inclusion=FALSE, warnOnly=TRUE)
checkNumeric(pi, 1, 4, integer=TRUE, warnOnly=TRUE)
checkNumeric(c(1, 1), 1, 4, warnOnly=TRUE)
checkNumeric(c(1, 1), 1, 4, unique=FALSE, warnOnly=TRUE)

##
## checkLogicalInteger
##
checkLogicalInteger(NULL, 3)
checkLogicalInteger(c(FALSE, TRUE), warnOnly=TRUE) 
checkLogicalInteger(1:2, 3) 
checkLogicalInteger(2, warnOnly=TRUE) 
checkLogicalInteger(c(2, 4), 3, warnOnly=TRUE)

##
## checkLogicalInteger names its calling function 
## rather than itself as the location of error detection
## if possible
##
tstFun <- function(x, length., warnOnly=FALSE){
   checkLogicalInteger(x, length., warnOnly) 
}
tstFun(NULL, 3)
tstFun(4, 3, warnOnly=TRUE)

tstFun2 <- function(x, length., warnOnly=FALSE){
   tstFun(x, length., warnOnly)
}
tstFun2(4, 3, warnOnly=TRUE)

}
\keyword{attribute} 
\keyword{utilities}

