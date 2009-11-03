\name{pca.fd}
\alias{pca.fd}
\title{
  Functional Principal Components Analysis
}
\description{
  Functional Principal components analysis aims to display types of
  variation across a sample of functions.  Principal components analysis
  is an exploratory data analysis that tends to be an early part of many
  projects.  These modes of variation are called $principal components$
  or $harmonics.$  This function computes these harmonics, the
  eigenvalues that indicate how important each mode of variation, and
  harmonic scores for individual functions. If the functions are
  multivariate, these harmonics are combined into a composite function
  that summarizes joint variation among the several functions that make
  up a multivariate functional observation.
}
\usage{
pca.fd(fdobj, nharm = 2, harmfdPar=fdPar(fdobj),
       centerfns = TRUE)
}
\arguments{
  \item{fdobj}{
    a functional data object.
  }
  \item{nharm}{
    the number of harmonics or principal components to compute.
  }
  \item{harmfdPar}{
    a functional parameter object that defines the
    harmonic or principal component functions to be estimated.
  }
  \item{centerfns}{
    a logical value:
    if TRUE, subtract the mean function from each function before
    computing principal components.
  }
}
\value{
  an object of class "pca.fd" with these named entries:

  \item{harmonics}{
    a functional data object for the harmonics or eigenfunctions
  }
  \item{values}{
    the complete set of eigenvalues
  }
  \item{scores}{
    s matrix of scores on the principal components or harmonics
  }
  \item{varprop}{
    a vector giving the proportion of variance explained
    by each eigenfunction
  }
  \item{meanfd}{
    a functional data object giving the mean function
  }
}
\seealso{
\code{\link{cca.fd}},
\code{\link{pda.fd}}
}
\examples{

#  carry out a PCA of temperature
#  penalize harmonic acceleration, use varimax rotation

daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)

harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
harmfdPar     <- fdPar(daybasis65, harmaccelLfd, lambda=1e5)
daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
                   daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd

daytemppcaobj <- pca.fd(daytempfd, nharm=4, harmfdPar)
daytemppcaVarmx <- varmx.pca.fd(daytemppcaobj)
#  plot harmonics
op <- par(mfrow=c(2,2))
plot.pca.fd(daytemppcaobj, cex.main=0.9)

plot.pca.fd(daytemppcaVarmx, cex.main=0.9)
par(op)

plot(daytemppcaobj$harmonics)
plot(daytemppcaVarmx$harmonics)
}
% docclass is function
\keyword{smooth}
