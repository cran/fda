\name{knots.fd}
\alias{knots.fd}
\alias{knots.fdSmooth}
\alias{knots.basisfd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract the knots from a function basis or data object}
\description{
  Extract either all or only the interior knots from an object of class
  \code{basisfd}, \code{fd}, or \code{fdSmooth}.
}
\usage{
\method{knots}{fd}(Fn, interior=TRUE, ...)
\method{knots}{fdSmooth}(Fn, interior=TRUE, ...)
\method{knots}{basisfd}(Fn, interior=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fn}{
    an object of class \code{basisfd} or containing such an object
  }
  \item{interior}{
    logical:

    if TRUE, Fn[["params"]] are returned.

    Else, nord <- norder(Fn); rng <- Fn[['rangeval']]; return
    c(rep(rng[1], nord), Fn[["params"]], rep(rng[2], nord))
  }
  \item{\dots}{ignored}
}
\details{
  The interior knots of a \code{bspline} basis are stored in the
  \code{params} component.  The remaining knots are in the
  \code{rangeval} component, with multiplicity norder(Fn).
}
\value{
  Numeric vector.  If 'interior' is TRUE, this is the \code{params}
  component of the \code{bspline} basis.  Otherwise, \code{params} is
  bracketed by rep(rangeval, norder(basisfd).
}
\author{Spencer Graves }
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{fd}},
  \code{\link{create.bspline.basis}},
}
% \examples{
% x <- 0:24
% y <- c(1.0,1.0,1.4,1.1,1.0,1.0,4.0,9.0,13.0,
       % 13.4,12.8,13.1,13.0,14.0,13.0,13.5,
       % 10.0,2.0,3.0,2.5,2.5,2.5,3.0,4.0,3.5)

% # knots.fdSmooth
% girlGrowthSm <- with(growth, smooth.basisPar(argvals=age, y=hgtf, lambda=0.1)$fd)

% girlKnots.fdSm  <- knots(girlGrowthSm)
% girlKnots.fdSmA <- knots(girlGrowthSm, interior=FALSE)
% stopifnot(all.equal(girlKnots.fdSm, girlKnots.fdSmA[5:33]))

% girlKnots.fd  <- knots(girlGrowthSm)
% girlKnots.fdA <- knots(girlGrowthSm, interior=FALSE)

% stopifnot(all.equal(girlKnots.fdSm,  girlKnots.fd))
% stopifnot(all.equal(girlKnots.fdSmA, girlKnots.fdA))

% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{optimize}
