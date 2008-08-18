\name{fdlabels}
\alias{fdlabels}
\title{
  Extract plot labels and names for replicates and variables 
}

\description{
  Extract plot labels and, if available, names for each replicate and
  variable 
}
\usage{
fdlabels(fdnames, nrep, nvar) 
}
\arguments{
  \item{fdnames}{
    a list of lenght 3 with xlabel, casenames, and ylabels.  
  }
  \item{nrep}{ integer number of cases or observations }
  \item{nvar}{integer number of variables }
}
\details{
  xlabel <- if(length(fdnames[[1]])>1) names(fdnames)[1] else 
  fdnames[[1]] 

  ylabel <- if(length(fdnames[[3]])>1) names(fdnames)[3] else
  fdnames[[3]] 

  casenames <- if(length(fdnames[[2]])== nrep)fdnames[[2]] else NULL

  varnames <- if(length(fdnames[[3]])==nvar)fdnames[[3]] else NULL 
}
\value{
  A list of xlabel, ylabel, casenames, and varnames 
}
\seealso{
  \code{\link{plot.fd}} 
}
\author{Jim Ramsay}
%\examples{}
\keyword{smooth }

