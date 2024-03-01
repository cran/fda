\name{cumfd}
\alias{cumfd}
\title{
  Compute a Cumulative Distribution Functional Data Object
}
\description{
  Function \code{smooth.morph()} maps a sorted set of variable values inside
  a closed interval into a set of equally-spaced probabilities in [0,1].
}
\usage{
cumfd(xrnd, xrng, nbreaks=7, nfine=101)
}
\arguments{
  \item{xrnd}{A vector of variable unsorted values.}
  \item{xrng}{A vector of length 2 containing the boundary values.}
  \item{nbreaks}{The number of knots to use to define object \code{WfdPar} in
  function \code{smooth.morph()}.}
  \item{nfine}{The number of equally spaced values spanning xrng.}
}
\details{
  Only the values of x within the interior of xrng are used 
  in order to avoid distortion due to boundary inflation or deflation.
}
\value{
  A named list of length 2 containing:
  \item{Wfdobj}{
    a functional data object defining function $W(x)$ that that
    optimizes the fit to the data of the monotone function that it
    defines.
  }
  \item{cdffine}{
    a vector of length nfine of an equally spaced mesh of values for
    the cumulative distribution function.
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
  \code{\link{smooth.morph}},
  \code{\link{landmarkreg}},
  \code{\link{register.fd}}
}
\examples{
#  see the use of smooth.morph in landmarkreg.R
xrnd <- rbeta(50, 2, 5)
xrng <- c(0,1)
hist(xrnd)
range(xrnd)
cdfd <- cumfd(xrnd, xrng)
}
