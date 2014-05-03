\name{monfn}
\alias{monfn}
\title{Evaluates a monotone function}
\description{Evaluates a monotone function}
\usage{
  monfn(argvals, Wfdobj, basislist=vector("list", JMAX), returnMatrix=FALSE)
}
\arguments{
  \item{argvals}{A numerical vector at which function and derivative are
                 evaluated.}
  \item{Wfdobj}{A functional data object.}
  \item{basislist}{A list containing values of basis functions.}
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
}
\value{
  A numerical vector or matrix containing the values the warping function h.
}
\details{
  This function evaluates a strictly monotone function of the form
  \deqn{h(x) = [D^{-1} exp(Wfdobj)](x),}
  where \eqn{D^{-1}} means taking the indefinite integral. The interval over
  which the integration takes places is defined in the basis object in Wfdobj.
}
\seealso{
  \code{\link{landmarkreg}}
}
\examples{

## basically this example resembles part of landmarkreg.R that uses monfn.R to
## estimate the warping function.

## Specify the curve subject to be registered
n=21
tbreaks = seq(0, 2*pi, len=n)
xval <- sin(tbreaks)
rangeval <- range(tbreaks)

## Establish a B-spline basis for the curve
wbasis <- create.bspline.basis(rangeval=rangeval, breaks=tbreaks)
Wfd0   <- fd(matrix(0,wbasis$nbasis,1),wbasis)
WfdPar <- fdPar(Wfd0, 1, 1e-4)
fdObj  <- smooth.basis(tbreaks, xval, WfdPar)$fd

## Set the mean landmark times. Note that the objective of the warping
## function is to transform the curve such that the landmarks of the curve
## occur at the designated mean landmark times.

## Specify the mean landmark times: tbreak[8]=2.2 and tbreaks[13]=3.76
meanmarks <- c(rangeval[1], tbreaks[8], tbreaks[13], rangeval[2])
## Specify landmark locations of the curve: tbreaks[6] and tbreaks[16]
cmarks <- c(rangeval[1], tbreaks[6], tbreaks[16], rangeval[2])

## Establish a B-basis object for the warping function
Wfd = smooth.morph(x=meanmarks, y=cmarks, WfdPar=WfdPar)$Wfdobj

## Estimate the warping function
h = monfn(tbreaks, Wfd)

## scale using a linear equation h such that h(0)=0 and h(END)=END
b <- (rangeval[2]-rangeval[1])/ (h[n]-h[1])
a <- rangeval[1] - b*h[1]
h <- a + b*h
plot(tbreaks, h, xlab="Time", ylab="Transformed time", type="b")
}

