\name{cycleplot.fd}
\alias{cycleplot.fd}
\title{
  Plot Cycles for a Periodic Bivariate Functional Data Object
}
\description{
A plotting function for data such as the knee-hip angles in
the gait data or temperature-precipitation curves for the
weather data.
}
\usage{
cycleplot.fd(fdobj, matplt=TRUE, nx=201, ...)
}
\arguments{
\item{fdobj}{
a bivariate functional data object to be plotted.
}
\item{matplt}{
if TRUE, all cycles are plotted simultaneously; otherwise
each cycle in turn is plotted.
}
\item{nx}{
the number of argument values in a fine mesh to be
plotted.  Increase the default number of 201 if
the curves have a lot of detail in them.
}
\item{\dots }{
additional plotting parameters such as axis labels and
etc. that are used in all plot functions.
}
}
\value{
None
}
\section{Side Effects}{
A plot of the cycles
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
  \code{\link{plot.fd}}, 
  \code{\link{plotfit.fd}}, 
}
% docclass is function
\keyword{smooth}
