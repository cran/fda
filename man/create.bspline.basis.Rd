\name{create.bspline.basis}
\alias{create.bspline.basis}
\title{Create a set of B-spline basis functions for constructing a functional data object. }
\description{
Functional data objects are constructed by specifying a set of basis functions and a set of coefficients defining a linear combination of these basis functions. The B-spline basis is a system that is usually used for non-periodic functions. It has the advantages of very fast computation and great flexibility. On the other hand, if the data are considered to be periodic, then the Fourier basis is usually preferred. 
}
\usage{
create.bspline.basis(rangeval=c(0, 1), nbasis=NULL, norder=4, breaks=NULL)
}
\arguments{
  \item{rangeval}{A vector of length 2 containing the initial and final values of the interval over which the functional data object can be evaluated. Default value c(0,1) }
  \item{nbasis}{An integer variable specifying the number of basis functions. Default value NULL.}
  \item{norder}{An integer specifying the order of b-splines, which is one higher than their degree. The default of 4 gives cubic splines. }
  \item{breaks}{A vector specifying the break points defining the b-spline. Also called knots, these are a strictly increasing sequence of junction points between piecewise polynomial segments. They must satisfy breaks[1] = rangeval[1] and breaks[nbreaks] = rangeval[2], where nbreaks is the length of breaks. There must be at least 3 values in breaks. There is a potential for inconsistency among arguments nbasis, norder, and breaks. It is resolved as follows: If breaks is supplied, nbreaks = length(breaks), and nbasis = nbreaks + norder - 2, no matter what value for nbasis is supplied. If breaks is not supplied, but nbasis is, nbreaks = nbasis - norder + 2, and if this turns out to be less than 3, an error message results. If neither breaks nor nbasis is supplied, nbreaks is set to 21. }
}
\details{
Spline functions are constructed by joining polynomials end-to-end at argument values called break points or knots . First, the interval is subdivided into a set of adjoining intervals separated the knots. Then a polynomial of order m and degree m-1 is defined for each interval. In order to make the resulting piece-wise polynomial smooth, two adjoining polynomials are constrained to have their values and all their derivatives up to order m-2 match at the point where they join. 

Consider as an illustration the very common case where the order is 4 for all polynomials, so that degree of each polynomials is 3. That is, the polynomials are cubic . Then at each break point or knot, the values of adjacent polynomials must match, and so also for their first and second derivatives. Only their third derivatives will differ at the point of junction. 

The number of degrees of freedom of a cubic spline function of this nature is calculated as follows. First, for the first interval, there are four degrees of freedom. Then, for each additional interval, the polynomial over that interval now has only one degree of freedom because of the requirement for matching values and derivatives. This works out to the following equation: The number of degrees of freedom is equal to the number of interior knots (that is, not counting the lower and upper limits) plus the order of the polynomials. 

B-splines are a set of special spline functions that can be used to construct any such piece-wise polynomial by computing the appropriate linear combination. They derive their computational convience from the fact that any B-spline basis function is nonzero over at most m adjacent intervals. The number of basis functions is given by the rule above for the number of degrees of freedom. 

The number of intervals controls the flexibility of the spline; the more knots, the more flexible the resulting spline will be. But the position of the knots also plays a role. Where do we position the knots? There is room for judgment here, but two considerations must be kept in mind: (1) you usually want at least one argument value between two adjacent knots, and (2) there should be more knots where the curve needs to have sharp curvatures such as a sharp peak or valley or an abrupt change of level, but only a few knots are required where the curve is changing very slowly. 

}
\value{
Returns a list with the basis.fd, class attribute with entries as above having names type, rangeval, nbasis, and params, respectively. The type member is bspline. 
}
\references{}
\note{}
\author{}
\seealso{
basis.fd. create.basis, create.fourier.basis 
}
\examples{}
\keyword{}
