\name{basisfd.product}
\alias{basisfd.product}
\alias{*.basisfd}
\title{
  Product of two basisfd objects
}
\description{
   pointwise multiplication method for basisfd class
}
\usage{
\method{*}{basisfd}(basisobj1, basisobj2)
}
\arguments{
\item{basisobj1, basisobj2}{
  objects of class basisfd
}
}
\details{
  TIMES for (two basis objects sets up a basis suitable for expanding
  the pointwise product of two functional data objects with these
  respective bases.  In the absence of a true product basis system in
  this code, the rules followed are inevitably a compromise:
   (1) if both bases are B-splines, the norder is the sum of the
       two orders - 1, and the breaks are the union of the
       two knot sequences, each knot multiplicity being the maximum
       of the multiplicities of the value in the two break sequences.
       That is, no knot in the product knot sequence will have a
       multiplicity greater than the multiplicities of this value
       in the two knot sequences.
       The rationale this rule is that order of differentiability
       of the product at eachy value will be controlled  by
       whichever knot sequence has the greater multiplicity.
       In the case where one of the splines is order 1, or a step
       function, the problem is dealt with by replacing the
       original knot values by multiple values at that location
       to give a discontinuous derivative.
   (2) if both bases are Fourier bases, AND the periods are the
       the same, the product is a Fourier basis with number of
       basis functions the sum of the two numbers of basis fns.
   (3) if only one of the bases is B-spline, the product basis
       is B-spline with the same knot sequence and order two
       higher.
   (4) in all other cases, the product is a B-spline basis with
       number of basis functions equal to the sum of the two
       numbers of bases and equally spaced knots.

}
\seealso{
  \code{\link{basisfd}}
}
\examples{
  f1 <- create.fourier.basis()
  f1.2 <- f1*f1
  \dontshow{stopifnot(}
  all.equal(f1.2, create.fourier.basis(nbasis=5))
  \dontshow{)}
}
% docclass is function
\keyword{smooth}

