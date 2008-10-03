\name{AmpPhaseDecomp}
\alias{AmpPhaseDecomp}
\title{
  Decomposition for Amplitude and Phase Variation
}
\description{
  Registration is the process of aligning peaks, valleys and other
  features in a sample of curves.  Once the registration has taken
  place, this function computes two mean squared error measures, one for
  amplitude variation, and the other for phase variation.  It also
  computes a squared multiple correlation index of the amount of
  variation in the unregistered functions is due to  phase.
}
\usage{
AmpPhaseDecomp(xfd, yfd, wfd)
}
\arguments{
  \item{xfd}{
    a functional data object containing the unregistered curves.
  }
  \item{yfd}{
    a functional data object containing the registered curves.
  }
  \item{wfd}{
    a functional data object defining the functions $w(t)$ that define
    the strictly monotone warping functions $h(t)$.  This is typically
    returned by the functions \code{landmarkreg} and
    \code{register.fd}.
  }
}
\value{
  a named list with the following components:

  \item{MS.amp}{
    the mean squared error for amplitude variation.
  }
  \item{MS.phas}{
    the mean squared error for phase variation.
  }
  \item{RSQR}{
    the squared correlation measure of the  proportion of the total
    variation that is due to phase variation.
  }
  \item{C}{
    a constant required for the decomposition.  Its value is one if the
    derivatives the warping functions are independent of the squared
    registered functions.
  }
}
\details{
  The decomposition can yield negative values for \code{MS.phas} if the
  registration does not improve the alignment of the curves, or if used
  to compare two registration processes based on different principles,
  such as  is the case for functions \code{landmarkreg} and
  \code{register.fd}.
}
\seealso{
  \code{\link{landmarkreg}},
  \code{\link{register.fd}},
  \code{\link{smooth.morph}}
}
\examples{
#See the analysis for the growth data in the examples.
}
% docclass is function
\keyword{smooth}
\keyword{models}

