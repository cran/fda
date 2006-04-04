\name{lines.fd}
\alias{lines.fd}
\title{
  Add Lines from Functional Data to a Plot
}
\description{
Lines defined by functional observations are added to an existing plot.
}
\usage{
lines.fd(x, Lfdobj=int2Lfd(0), ...)
}
\arguments{
\item{x}{
a univariate functional data object.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
If present, the derivative or the value of applying the operator is
evaluated rather than the functions themselves.
}
\item{\dots }{
additional arguments such as axis titles and so forth that can be used
in all plotting programs.
}
}
\section{Side Effects}{
Lines added to an existing plot.
}
\seealso{
plot.fd, plotfit.fd
}
% docclass is function
\keyword{smooth}
