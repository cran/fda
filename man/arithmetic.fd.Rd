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
"+.fd"(e1, e2)
"-.fd"(e1, e2)
"*.fd"(e1, e2)
plus.fd(e1, e2, basisobj=basisobj1*basisobj2)
minus.fd(e1, e2, basisobj=basisobj1*basisobj2)
times.fd(e1, e2, basisobj=basisobj1*basisobj2)
}
\arguments{
  \item{e1, e2}{
    object of class 'fd' or a numeric vector.  Note that 'e1+e2' will
    dispatch to plus.fd(e1, e2) only of e1 has class 'fd'.  Similarly,
    'e1-e2' or 'e1*e2' will dispatch to minus.fd(e1, e2) or time.fd(e1,
    e2), respetively, only if e1 is of class 'fd'.
  }
  \item{basisobj}{
    reference basis
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
\code{\link{basisfd.product}}
}
\examples{
##
## add a parabola to itself
##
bspl4 <- create.bspline.basis(nbasis=4)
parab4.5 <- fd(c(3, -1, -1, 3)/3, bspl4)
\dontshow{stopifnot(}
all.equal(coef(parab4.5+parab4.5), matrix(c(6, -2, -2, 6)/3, 4))
\dontshow{)}
\dontshow{stopifnot(}
all.equal(coef(parab4.5-parab4.5), matrix(rep(0, 4), 4))
\dontshow{)}

##
## Same example with interior knots at 1/3 and 1/2
##
bspl5.3 <- create.bspline.basis(breaks=c(0, 1/3, 1))
plot(bspl5.3)
x. <- seq(0, 1, .1)
para4.5.3 <- smooth.basis(x., 4*(x.-0.5)^2, fdParobj=bspl5.3)[['fd']]
plot(para4.5.3)

bspl5.2 <- create.bspline.basis(breaks=c(0, 1/2, 1))
plot(bspl5.2)
para4.5.2 <- smooth.basis(x., 4*(x.-0.5)^2, fdParobj=bspl5.2)[['fd']]
plot(para4.5.2)

str(para4.5.3+para4.5.2)
\dontshow{stopifnot(}
all.equal(coef(para4.5.3-para4.5.2), matrix(0, 9, 1))
\dontshow{)}

str(para4.5.3*para4.5.2)
# interior knots of the sum
# = union(interior knots of the summands);
# ditto for difference and product.
plot(para4.5.3*para4.5.2)

##
## fd+numeric
##
\dontshow{stopifnot(}
all.equal(coef(parab4.5+1), matrix(c(6, 2, 2, 6)/3, 4))
\dontshow{)}

\dontshow{stopifnot(}
all.equal(1+parab4.5, parab4.5+1)
\dontshow{)}

##
## fd-numeric
##
\dontshow{stopifnot(}
all.equal(coef(-parab4.5), matrix(c(-3, 1, 1, -3)/3, 4))
\dontshow{)}


plot(parab4.5-1)

plot(1-parab4.5)

}
% docclass is function
\keyword{smooth}

