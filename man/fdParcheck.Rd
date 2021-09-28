\name{fdParcheck}
\alias{fdParcheck}
\title{
Convert \code{fd} or \code{basisfd} Objects to \code{fdPar} Objects
}
\description{
If the input is an \code{fd} object, default parameters are added to convert
it to an \code{fdPar} object.  If a \code{basisfd} object,  it is first
converted to an \code{fd} object with a square coefficient matrix of
zeros, and then converted to a \code{fdPar} object.
}
\usage{
fdParcheck(fdPar)
}
\arguments{
  \item{fdPar}{A \code{basisfd}, \code{fd} or a \code{fdPar} object.}
}
\value{
Returns a \code{fdPar} object.
}
