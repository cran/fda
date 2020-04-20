\name{phaseplanePlot}
\alias{phaseplanePlot}
\title{
  Phase-plane plot
}
\description{
  Plot acceleration (or Ldfobj2) vs. velocity (or Lfdobj1) of a function
  data object.
}
\usage{
phaseplanePlot(evalarg, fdobj, Lfdobj1=1, Lfdobj2=2,
        lty=c("longdash", "solid"),
      labels=list(evalarg=seq(evalarg[1], max(evalarg), length=13),
             labels=fda::monthLetters),
      abline=list(h=0, v=0, lty=2), xlab="Velocity",
      ylab="Acceleration", returnMatrix=FALSE, ...)
}
\arguments{
  \item{evalarg}{
    a vector of argument values at which the functional data object is
    to be evaluated.

    Defaults to a sequence of 181 points in the range
    specified by fdobj[["basis"]][["rangeval"]].

    If(length(evalarg) == 1)it is replaced by seq(evalarg[1],
    evalarg[1]+1, length=181).

    If(length(evalarg) == 2)it is replaced by seq(evalarg[1],
    evalarg[2], length=181).
  }
  \item{fdobj}{
    a functional data object to be evaluated.
  }
  \item{Lfdobj1}{
    either a nonnegative integer or a linear differential operator
    object.  The points plotted on the horizontal axis are
    eval.fd(evalarg, fdobj, Lfdobj1).  By default, this is the
    velocity.
  }
  \item{Lfdobj2}{
    either a nonnegative integer or a linear differential operator
    object.  The points plotted on the vertical axis are
    eval.fd(evalarg, fdobj, Lfdobj2).  By default, this is the
    acceleration.
  }
  \item{lty}{
    line types for the first and second halves of the plot.
  }
  \item{labels}{
    a list of length two:

    evalarg = a numeric vector of 'evalarg' values to be labeled.

    labels = a character vector of labels, replicated to the same length
    as labels[["evalarg"]] in case it's not of the same length.
  }
  \item{abline}{
    arguments to a call to abline.
  }
  \item{xlab}{
    x axis label
  }
  \item{ylab}{
    y axis label
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
  \item{...}{
    optional arguments passed to plot.
  }

}
\value{
  Invisibly returns a matrix with two columns containg the points
  plotted.
}
\seealso{
  \code{\link{plot}},
  \code{\link{eval.fd}}
  \code{\link{plot.fd}}
  \code{\link{nondurables}}
}
\examples{
goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=161, norder=8)
LfdobjNonDur <- int2Lfd(4)
argvals = seq(1919,2000,len=length(nondurables))
logNondurSm <- smooth.basisPar(argvals,
                y=log10(nondurables), fdobj=goodsbasis,
                Lfdobj=LfdobjNonDur, lambda=1e-11)
phaseplanePlot(1964, logNondurSm$fd)
}
\keyword{smooth}
\keyword{hplot}
