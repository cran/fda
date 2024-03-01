\name{plot.pda.fd}
\alias{plot.pda.fd}
\title{
  Plot Principle Differential Analysis Components
}
\description{
  Plots the results of pda.fd, allows the user to group coefficient functions
  by variable, equation, derivative or combination of them.
}
\usage{
%plot.pda(pdaList,whichdim=1,npts=501,...)
\method{plot}{pda.fd}(x,whichdim=1,npts=501,...)
}
\arguments{
  \item{x}{
    an object of class \code{pda.fd}.
  }
  \item{whichdim}{
    which dimension to use as grouping variables
    \describe{
      \item{1}{ coefficients of each variable differential equation}
      \item{2}{ coefficient functions for each equation}
      \item{3}{ coefficients of derivatives of each variable}
    }
    \code{whichdim} should be an ordered vector of length between 1 and
    3.
  }
  \item{npts}{
    number of points to use for plotting.
  }
  \item{\dots}{
    other arguments for 'plot'.
  }
}
\details{
  Produces one plot for each coefficient function in a principle differential
  analysis.
}
\value{
  invisible(NULL)
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
  \code{\link{pda.fd}}
  \code{\link{eigen.pda}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
#  A pda analysis of the handwriting data

# reduce the size to reduce the compute time for the example
ni <- 281
indx <- seq(1, 1401, length=ni)
fdaarray <- handwrit[indx,,]
fdatime  <- seq(0, 2.3, len=ni)

#  basis for coordinates

fdarange <- c(0, 2.3)
breaks   <- seq(0,2.3,length.out=116)
norder   <- 6
fdabasis <- create.bspline.basis(fdarange,norder=norder,breaks=breaks)

#  parameter object for coordinates

fdafd0 <- fd(matrix(0,fdabasis$nbasis,1), fdabasis)
fdaPar <- fdPar(fdafd0,int2Lfd(4),1e-8)

#  coordinate functions and a list tontaining them

Xfd <- smooth.basis(fdatime, fdaarray[,,1], fdaPar)$fd
Yfd <- smooth.basis(fdatime, fdaarray[,,2], fdaPar)$fd

xfdlist <- list(Xfd, Yfd)

#  basis and parameter object for weight functions

fdabasis2 <- create.bspline.basis(fdarange,norder=norder,nbasis=31)
fdafd0    <- fd(matrix(0,fdabasis2$nbasis,1), fdabasis2)
pdaPar    <- fdPar(fdafd0,1,1e-8)

pdaParlist <- list(pdaPar, pdaPar)

bwtlist <- list( list(pdaParlist,pdaParlist), list(pdaParlist,pdaParlist) )

#  do the second order pda

pdaList <- pda.fd(xfdlist, bwtlist)

# plot the results

plot(pdaList,whichdim=1)
plot(pdaList,whichdim=2)
plot(pdaList,whichdim=3)

plot(pdaList,whichdim=c(1,2))
plot(pdaList,whichdim=c(1,3))
plot(pdaList,whichdim=c(2,3))

plot(pdaList,whichdim=1:3)
par(oldpar)
}

\keyword{smooth}
