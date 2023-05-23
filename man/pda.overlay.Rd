\name{pda.overlay}
\alias{pda.overlay}
\title{
  Stability Analysis for Principle Differential Analysis
}
\description{
  Overlays the results of a univariate, second-order principal differential
  analysis on a bifurcation diagram to demonstrate stability.
}
\usage{
pda.overlay(pdaList,nfine=501,ncoarse=11,...)
}
\arguments{
  \item{pdaList}{
    a list object returned by \code{pda.fd}.
  }
  \item{nfine}{
    number of plotting points to use.
  }
  \item{ncoarse}{
    number of time markers to place along the plotted curve.
  }
  \item{\dots}{
    other arguments for 'plot'.
  }
}
\details{
  Overlays a bivariate plot of the functional parameters in a univariate
  second-order principal differential analysis on a bifurcation diagram.
}
\value{
  None.
}
\seealso{
  \code{\link{pda.fd}}
  \code{\link{plot.pda.fd}}
  \code{\link{eigen.pda}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
# This example looks at a principal differential analysis of the lip data
# in Ramsay and Silverman (2005).

# First smooth the data

lipfd <- smooth.basisPar(liptime, lip, 6, Lfdobj=int2Lfd(4),
                         lambda=1e-12)$fd
names(lipfd$fdnames) <- c("time(seconds)", "replications", "mm")

# Now we'll set up functional parameter objects for the beta coefficients.

lipbasis <- lipfd$basis
lipfd0   <- fd(matrix(0,lipbasis$nbasis,1),lipbasis)
lipfdPar <- fdPar(lipfd0,2,0)
bwtlist  <- list(lipfdPar,lipfdPar)
xfdlist  <- list(lipfd)

# Call pda

pdaList <- pda.fd(xfdlist, bwtlist)

# And plot the overlay

pda.overlay(pdaList,lwd=2,cex.lab=1.5,cex.axis=1.5)
par(oldpar)
}
\keyword{smooth}
