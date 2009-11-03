\name{lmWinsor12}
\alias{lmWinsor1}
\alias{lmWinsor2}
\alias{constraintCheck}
\title{
  Support functions for lmWinsor
}
\description{
  internal functions used by lmWinsor, not intended for direct use and
  may be changed in the future without notice.
}
\usage{
constraintCheck(Amat, bvec, eps=sqrt(.Machine$double.eps) )
lmWinsor1(formula, data, lower=NULL, upper=NULL, trim=0,
        quantileType=7, subset, weights=NULL, na.action,
        model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
        singular.ok = TRUE, contrasts = NULL, offset=NULL,
        method=c('QP', 'clip'), eps=sqrt(.Machine$double.eps),
        trace=ceiling(sqrt(nrow(data))), ...)
lmWinsor2(formula, data, lower=NULL, upper=NULL, trim=0,
        quantileType=7, subset, weights=NULL, na.action,
        model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
        singular.ok = TRUE, contrasts = NULL, offset=NULL,
        method=c('QP', 'clip'), eps=sqrt(.Machine$double.eps),
        trace=ceiling(sqrt(nrow(data))), ...)
}
\arguments{
  \item{formula, data, lower, upper, trim, quantileType, subset,
    weights, na.action, model, x, y, qr, singular.ok, contrasts,
    offset, method, eps, trace, \dots}{as for lmWinsor}
  \item{Amat, bvec}{t(Amat) and bvec are arguments to \code{solve.QP}.}
}
%\details{}
%\value{}
\author{ Spencer Graves }
\seealso{
  \code{\link{lmWinsor}}
  \code{\link[quadprog]{solve.QP}}
}
%\examples{}
\keyword{ models }
