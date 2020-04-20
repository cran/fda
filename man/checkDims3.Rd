\name{checkDims3}
\alias{checkDim3}
\alias{checkDims3}
\title{
  Compare dimensions and dimnames of arrays
}
\description{
  Compare selected dimensions and dimnames of arrays, coercing objects
  to 3-dimensional arrays and either give an error or force matching.
}
\usage{
checkDim3(x, y=NULL, xdim=1, ydim=1, defaultNames='x',
         subset=c('xiny', 'yinx', 'neither'),
         xName=substring(deparse(substitute(x)), 1, 33),
         yName=substring(deparse(substitute(y)), 1, 33) )
checkDims3(x, y=NULL, xdim=2:3, ydim=2:3, defaultNames='x',
         subset=c('xiny', 'yinx', 'neither'),
         xName=substring(deparse(substitute(x)), 1, 33),
         yName=substring(deparse(substitute(y)), 1, 33) )
}
\arguments{
  \item{x, y}{
    arrays to be compared.  If \code{y} is missing, \code{x} is used.

    Currently, both \code{x} and \code{y} can have at most 3
    dimensions.  If either has more, an error will be thrown.  If either
    has fewer, it will be expanded to 3 dimensions using
    \code{as.array3}.
  }
  \item{xdim, ydim}{
    For \code{checkDim3}, these are positive integers indicating which
    dimension of \code{x} will be compared with which dimension of
    \code{y}.

    For \code{checkDims3}, these are positive integer vectors of the same
    length, passed one at a time to \code{checkDim3}.  The default here
    is to force matching dimensions for \code{\link{plotfit.fd}}.
  }
  \item{defaultNames}{
    Either NULL, FALSE or a character string or vector or list.  If
    NULL, no checking is done of dimnames.  If FALSE, an error is thrown
    unless the corresponding dimensions of \code{x} and \code{y} match
    exactly.

    If it is a character string, vector, or list, it is used as the
    default names if neither \code{x} nor \code{y} have dimenames for
    the compared dimensions.  If it is a character vector that is too
    short, it is extended to the required length using
    paste(defaultNames, 1:ni), where \code{ni} = the required length.

    If it is a list, it should have length (length(xdim)+1).  Each
    component must be either a character vector or NULL.  If neither
    \code{x} nor \code{y} have dimenames for the first compared
    dimensions, defaultNames[[1]] will be used instead unless it is
    NULL, in which case the last component of defaultNames will be
    used.  If it is null, an error is thrown.
  }
  \item{subset}{
    If 'xiny', and any(dim(y)[ydim] < dim(x)[xdim]), an error is
    thrown.  Else if any(dim(y)[ydim] > dim(x)[xdim]) the larger is
    reduced to match the smaller.  If 'yinx', this procedure is
    reversed.

    If 'neither', any dimension mismatch generates an error.
  }

  \item{xName, yName}{
    names of the arguments \code{x} and \code{y}, used only to in error
    messages.
  }
}
\details{
  For \code{checkDims3}, confirm that \code{xdim} and \code{ydim} have
  the same length, and call \code{checkDim3} for each pair.

  For \code{checkDim3}, proceed as follows:

  1.  if((xdim>3) | (ydim>3)) throw an error.

  2.  ixperm <- list(1:3, c(2, 1, 3), c(3, 2, 1))[xdim];
  iyperm <- list(1:3, c(2, 1, 3), c(3, 2, 1))[ydim];

  3.  x3 <- aperm(as.array3(x), ixperm);
  y3 <- aperm(as.array3(y), iyperm)

  4.  xNames <- dimnames(x3);  yNames <- dimnames(y3)

  5.  Check subset.  For example, for subset='xiny', use the following:
  \code{
    if(is.null(xNames)){
      if(dim(x3)[1]>dim(y3)[1]) stop
      else y. <- y3[1:dim(x3)[1],,]
      dimnames(x) <- list(yNames[[1]], NULL, NULL)
    }
    else {
      if(is.null(xNames[[1]])){
	if(dim(x3)[1]>dim(y3)[1]) stop
	else y. <- y3[1:dim(x3)[1],,]
	dimnames(x3)[[1]] <- yNames[[1]]
      }
      else {
	if(any(!is.element(xNames[[1]], yNames[[1]])))stop
	else y. <- y3[xNames[[1]],,]
      }
    }
  }

  6.  return(list(x=aperm(x3, ixperm), y=aperm(y., iyperm)))
}
\value{
  a list with components \code{x} and \code{y}.
}
%\references{ }
\author{
  Spencer Graves
}
\seealso{
  \code{\link{as.array3}}
  \code{\link{plotfit.fd}}
}
\examples{
# Select the first two rows of y
stopifnot(all.equal(
checkDim3(1:2, 3:5),
list(x=array(1:2, c(2,1,1), list(c('x1','x2'), NULL, NULL)),
     y=array(3:4, c(2,1,1), list(c('x1','x2'), NULL, NULL)) )
))

# Select the first two rows of a matrix y
stopifnot(all.equal(
checkDim3(1:2, matrix(3:8, 3)),
list(x=array(1:2,         c(2,1,1), list(c('x1','x2'), NULL, NULL)),
     y=array(c(3:4, 6:7), c(2,2,1), list(c('x1','x2'), NULL, NULL)) )
))

# Select the first column of y
stopifnot(all.equal(
checkDim3(1:2, matrix(3:8, 3), 2, 2),
list(x=array(1:2,         c(2,1,1), list(NULL, 'x', NULL)),
     y=array(3:5, c(3,1,1), list(NULL, 'x', NULL)) )
))

# Select the first two rows and the first column of y
stopifnot(all.equal(
checkDims3(1:2, matrix(3:8, 3), 1:2, 1:2),
list(x=array(1:2, c(2,1,1), list(c('x1','x2'), 'x', NULL)),
     y=array(3:4, c(2,1,1), list(c('x1','x2'), 'x', NULL)) )
))

# Select the first 2 rows of y
x1 <- matrix(1:4, 2, dimnames=list(NULL, LETTERS[2:3]))
x1a <- x1. <- as.array3(x1)
dimnames(x1a)[[1]] <- c('x1', 'x2')
y1 <- matrix(11:19, 3, dimnames=list(NULL, LETTERS[1:3]))
y1a <- y1. <- as.array3(y1)
dimnames(y1a)[[1]] <- c('x1', 'x2', 'x3')

stopifnot(all.equal(
checkDim3(x1, y1),
list(x=x1a, y=y1a[1:2, , , drop=FALSE])
))

# Select columns 2 & 3 of y
stopifnot(all.equal(
checkDim3(x1, y1, 2, 2),
list(x=x1., y=y1.[, 2:3, , drop=FALSE ])
))

# Select the first 2 rows and  columns 2 & 3 of y
stopifnot(all.equal(
checkDims3(x1, y1, 1:2, 1:2),
list(x=x1a, y=y1a[1:2, 2:3, , drop=FALSE ])
))

# y = columns 2 and 3 of x
x23 <- matrix(1:6, 2, dimnames=list(letters[2:3], letters[1:3]))
x23. <- as.array3(x23)
stopifnot(all.equal(
checkDim3(x23, xdim=1, ydim=2),
list(x=x23., y=x23.[, 2:3,, drop=FALSE ])
))

# Transfer dimnames from y to x
x4a <- x4 <- matrix(1:4, 2)
y4 <- matrix(5:8, 2, dimnames=list(letters[1:2], letters[3:4]))
dimnames(x4a) <- dimnames(t(y4))
stopifnot(all.equal(
checkDims3(x4, y4, 1:2, 2:1),
list(x=as.array3(x4a), y=as.array3(y4))
))

# as used in plotfit.fd
daybasis65 <- create.fourier.basis(c(0, 365), 65)

daytempfd <- with(CanadianWeather, smooth.basis(
       day.5, dailyAv[,,"Temperature.C"], 
       daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd )

defaultNms <- with(daytempfd, c(fdnames[2], fdnames[3], x='x'))
subset <- checkDims3(CanadianWeather$dailyAv[, , "Temperature.C"],
               daytempfd$coef, defaultNames=defaultNms)
# Problem:  dimnames(...)[[3]] = '1'
# Fix:
subset3 <- checkDims3(
        CanadianWeather$dailyAv[, , "Temperature.C", drop=FALSE],
               daytempfd$coef, defaultNames=defaultNms)
}
\keyword{utilities}

