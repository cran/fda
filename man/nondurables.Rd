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

\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York. 

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York, ch. 3. 
}
\examples{
oldpar <- par(no.readonly=TRUE)
plot(nondurables, log="y")
par(oldpar)
} 
\keyword{datasets}
