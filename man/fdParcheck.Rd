\name{fdParcheck}
\alias{fdParcheck}
\title{
Convert \code{fd} or \code{basisfd} Objects to \code{fdPar} Objects
}
\description{
If the input is an \code{fd} object, default \code{fdPar} parameters are 
added to convert the object to an \code{fdPar} object.  
If a \code{basisfd} object,  it is first converted to an \code{fd} object 
with a \code{nbasis} by \code{ncurve} coefficient matrix of zeros, 
and then converted to a \code{fdPar} object.
}
\usage{
fdParcheck(fdParobj, ncurve=NULL)
}
\arguments{
  \item{fdParobj}{A \code{basisfd}, \code{fd} or a \code{fdPar} object.}
  \item{ncurve}{The number of curves or functions for which the basis
  object is to be used in their representation.  If missing, an error
  message is displayed and computation terminated.}
}
\details{
  The previous version of this function did not include the argument
  \code{ncurve}, and set up an \code{fd} object with a square
  matrix of zeros.  This could cause an error in many of the 
  functions that called it.
}
\value{
Returns a \code{fdPar} object.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
