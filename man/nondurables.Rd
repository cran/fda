\name{nondurables}
\alias{nondurables}
\title{
  Nondurable goods index
}
\description{
  US nondurable goods index time series, January 1919 to January 2000.  
}
\format{
    An object of class 'ts'.  
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\examples{
oldpar <- par(no.readonly=TRUE)
plot(nondurables, log="y")
par(oldpar)
} 
\keyword{datasets}
