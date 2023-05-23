\name{arithmetic.fd}
\alias{arithmetic.fd}
\alias{+.fd}
\alias{plus.fd}
\alias{-.fd}
\alias{minus.fd}
\alias{*.fd}
\alias{times.fd}
\title{
  Arithmetic on functional data ('fd') objects
}
\description{
  Arithmetic on functional data objects
}
\usage{
\method{+}{fd}(e1, e2)
\method{-}{fd}(e1, e2)
\method{*}{fd}(e1, e2)
plus.fd(e1, e2, basisobj=NULL)
minus.fd(e1, e2, basisobj=NULL)
times.fd(e1, e2, basisobj=NULL)
}
\arguments{
  \item{e1, e2}{
    object of class 'fd' or a numeric vector.  Note that 'e1+e2' will
    dispatch to plus.fd(e1, e2) only if e1 has class 'fd'.  Similarly,
    'e1-e2' or 'e1*e2' will dispatch to minus.fd(e1, e2) or time.fd(e1,
    e2), respectively, only if e1 is of class 'fd'.
  }
  \item{basisobj}{
    reference basis;  defaults to e1[['basis']] * e2[['basis']];
    ignored for \code{plus.fd} and \code{minus.fd}.
  }
}

\value{

  A function data object corresponding to the pointwise sum, difference
  or product of e1 and e2.

  If both arguments are functional data objects, the bases are the same,
  and the coefficient matrices are the same dims, the indicated
  operation is applied to the coefficient matrices of the two objects.
  In other words, e1+e2 is obtained for this case by adding the
  coefficient matrices from e1 and e2.

  If e1 or e2 is a numeric scalar,  that scalar is applied to the
  coefficient matrix of the functional data object.

  If either e1 or e2 is a numeric vector, it must be the same length as
  the number of replicated functional observations in the other
  argument.

  When both arguments are functional data objects, they need not have
  the same bases.  However, if they don't have the same number of
  replicates, then one of them must have a single replicate.  In the
  second case, the singleton function is replicated to match the number
  of replicates of the other function. In either case, they must have
  the same number of functions. When both arguments are functional data
  objects, and the bases are not the same, the basis used for the sum is
  constructed to be of higher dimension than the basis for either factor
  according to rules described in function TIMES for two basis objects.
}

\seealso{
  \code{\link{basisfd}},
  \code{\link{basisfd.product}},
  \code{\link{exponentiate.fd}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
##
## add a parabola to itself
##
bspl4 <- create.bspline.basis(nbasis=4)
parab4.5 <- fd(c(3, -1, -1, 3)/3, bspl4)

coef2 <- matrix(c(6, -2, -2, 6)/3, 4)
dimnames(coef2) <- list(NULL, 'reps 1')
\dontshow{stopifnot(}
all.equal(coef(parab4.5+parab4.5), coef2)
\dontshow{)}

##
## Same example with interior knots at 1/3 and 1/2
##
bspl5.3 <- create.bspline.basis(breaks=c(0, 1/3, 1))
plot(bspl5.3)
x. <- seq(0, 1, .1)
para4.5.3 <- smooth.basis(x., 4*(x.-0.5)^2, fdParobj=bspl5.3)$fd
plot(para4.5.3)

bspl5.2 <- create.bspline.basis(breaks=c(0, 1/2, 1))
plot(bspl5.2)
para4.5.2 <- smooth.basis(x., 4*(x.-0.5)^2, fdParobj=bspl5.2)$fd
plot(para4.5.2)

#str(para4.5.3+para4.5.2)

coef2. <- matrix(0, 9, 1)
dimnames(coef2.) <- list(NULL, 'rep1')
\dontshow{stopifnot(}
all.equal(coef(para4.5.3-para4.5.2), coef2.)
\dontshow{)}

##
## product
##
quart <- para4.5.3*para4.5.2

% # interior knots of the sum
% # = union(interior knots of the summands);
% # ditto for difference and product.
% \dontshow{stopifnot(}
% all.equal(knots.fd(quart), c(knots.fd(para4.5.3), knots.fd(para4.5.2)))
% \dontshow{)}

# norder(quart) = norder(para4.5.2)+norder(para4.5.3)-1 = 7
\dontshow{stopifnot(}
norder(quart) == (norder(para4.5.2)+norder(para4.5.3)-1)
\dontshow{)}

# para4.5.2 with knot at 0.5 and para4.5.3 with knot at 1/3
# both have (2 end points + 1 interior knot) + norder-2
#     = 5 basis functions
# quart has (2 end points + 2 interior knots)+norder-2
#     = 9 basis functions
# coefficients look strange because the knots are at
# (1/3, 1/2) and not symmetrical

\dontshow{stopifnot(}
all.equal(as.numeric(coef(quart)),
0.1*c(90, 50, 14, -10, 6, -2, -2, 30, 90)/9)
\dontshow{)}

plot(para4.5.3*para4.5.2) # quartic, not parabolic ...

##
## product with Fourier bases
##
f3 <- fd(c(0,0,1), create.fourier.basis())
f3^2 # number of basis functions = 7?

##
## fd+numeric
##
coef1 <- matrix(c(6, 2, 2, 6)/3, 4)
dimnames(coef1) <- list(NULL, 'reps 1')
\dontshow{stopifnot(}
all.equal(coef(parab4.5+1), coef1)
\dontshow{)}

\dontshow{stopifnot(}
all.equal(1+parab4.5, parab4.5+1)
\dontshow{)}

##
## fd-numeric
##
coefneg <- matrix(c(-3, 1, 1, -3)/3, 4)
dimnames(coefneg) <- list(NULL, 'reps 1')
\dontshow{stopifnot(}
all.equal(coef(-parab4.5), coefneg)
\dontshow{)}

plot(parab4.5-1)

plot(1-parab4.5)

par(oldpar)

}
\keyword{smooth}

