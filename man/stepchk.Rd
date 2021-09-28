\name{stepchk}
\alias{stepchk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Check a step size value for being within parameter limits.
}
\description{
Parameter values are check for a proposed stepsize to ensure that they
remain within limits.  If not, the stepsize is adjusted to keep within
the limits.}
\usage{
stepchk(oldstep, cvec, deltac, limwrd, ind, climit = 50 * c(-rep(1, ncvec), 
        rep(1, ncvec)), active = 1:ncvec, dbgwrd)
}
\arguments{
  \item{oldstep}{
A positive or zero step size for a line search}
  \item{cvec}{A numeric or matrix parameter vector.
}
  \item{deltac}{
A vector defining a research direction of the same length as \code{cvec}.}
  \item{limwrd}{
  A logical vector of length 2 indicating whether a lower and/or upper boundary is
  rached by a step.}
  \item{ind}{
  An integer that is normally 0 but is 1 if an upper boundary has been reached 
  twice.}
  \item{climit}{
  A two-column matrix with first dimension that of \code{cvec} indicating the lower   and upper boundaries for parameter vector values.}
  \item{active}{
  A logical vector with length that of \code{cvec} indicating which parmeter values   are fixed (FALSE) or free to be changed (TRUE).}
  \item{dbgwrd}{An integer controlling the amount of information to display for
  each step.}
}
\value{
A nonnegative step size for a linear search.
}