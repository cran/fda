\name{eigen.pda}
\alias{eigen.pda}
\title{
  Stability Analysis for Principle Differential Analysis
}
\description{
  Performs a stability analysis of the result of \code{pda.fd}, returning
  the real and imaginary parts of the eigenfunctions associated with the
  linear differential operator.
}
\usage{
eigen.pda(pdaList,plotresult=TRUE,npts=501,...)
}
\arguments{
  \item{pdaList}{
    a list object returned by \code{pda.fd}.
  }
  \item{plotresult}{
    should the result be plotted? Default is TRUE
  }
  \item{npts}{
    number of points to use for plotting.
  }
  \item{\dots}{
    other arguments for 'plot'.
  }
}
\details{
  Conducts an eigen decomposition of the linear differential equation implied
  by the result of \code{pda.fd}. Imaginary eigenvalues indicate instantaneous
  oscillatory behavior. Positive real eigenvalues indicate exponential increase,
  negative real eigenvalues correspond to exponential decay. If the principle
  differential analysis also included the estimation of a forcing function, the
  limitting stable points are also tracked.
}
\value{
  Returns a list with elements
  \item{argvals}{The evaluation points of the coefficient functions.}
  \item{eigvals}{The corresponding eigenvalues at each time.}
  \item{limvals}{The stable points of the system at each time.}
}
\seealso{
  \code{\link{pda.fd}}
  \code{\link{plot.pda.fd}}
  \code{\link{pda.overlay}}
}
\examples{

#  A pda analysis of the handwriting data

# reduce the size to reduce the compute time for the example
ni <- 281
indx <- seq(1, 1401, length=ni)
fdaarray = handwrit[indx,,]
fdatime  <- seq(0, 2.3, len=ni)

#  basis for coordinates

fdarange <- c(0, 2.3)
breaks = seq(0,2.3,length.out=116)
norder = 6
fdabasis = create.bspline.basis(fdarange,norder=norder,breaks=breaks)

#  parameter object for coordinates

fdaPar = fdPar(fdabasis,int2Lfd(4),1e-8)

#  coordinate functions and a list tontaining them

Xfd = smooth.basis(fdatime, fdaarray[,,1], fdaPar)$fd
Yfd = smooth.basis(fdatime, fdaarray[,,2], fdaPar)$fd

xfdlist = list(Xfd, Yfd)

#  basis and parameter object for weight functions

fdabasis2 = create.bspline.basis(fdarange,norder=norder,nbasis=31)
fdafd2    = fd(matrix(0,31,2),fdabasis2)
pdaPar    = fdPar(fdafd2,1,1e-8)

pdaParlist = list(pdaPar, pdaPar)

bwtlist = list( list(pdaParlist,pdaParlist), list(pdaParlist,pdaParlist) )

#  do the second order pda

pdaList = pda.fd(xfdlist, bwtlist)

# plot the results

eigres = eigen.pda(pdaList)
}
\keyword{smooth}
