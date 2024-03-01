\name{ppBspline}
\alias{ppBspline}
\title{Convert a B-spline function to piece-wise polynomial form}
\description{
  The B-spline basis functions of order \code{k = length(t) - 1}
  defined by the knot sequence in argument \code{t} each consist of polynomial
  segments with the same order joined end-to-end over the successive gaps in the
  knot sequence.  This function computes the \code{k} coefficients of these polynomial
  segments in the rows of the output matrix \code{coeff}, with each row corresponding
  to a B-spline basis function that is positive over the interval spanned by the
  values in \code{t}. The elements of the output vector \code{index} indicate where
  in the sequence \code{t} we find the knots.  Note that we assume
  \code{t[1] < t[k+1]}, i.e. \code{t} is not a sequence of the same knot.
}
\usage{
  ppBspline(t)
}
\arguments{
  \item{t}{
    numeric vector = knot sequence of length norder+1 where norder =
    the order of the B-spline.  The knot sequence must contain at least one gap.
  }
}
\value{
  a list object containing components
  \item{Coeff}{
    a matrix with rows corresponding to B-spline basis functions positive
    over the interval spanned by \code{t} and columns corresponding to the
    terms \code{1, x, x^2, ...} in the polynomial representation.
  }
  \item{index}{
    indices indicating where in the sequence \code{t} the knots are to be found
  }
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
    \code{\link{bsplineS}}
}
\examples{
  ppBspline(1:5)
}
\keyword{smooth}
