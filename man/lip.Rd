\name{lip}
\docType{data}
\alias{lip}
\alias{lipmarks}
\alias{liptime}
\title{
  Lip motion
}
\description{
  51 measurements of the position of the lower lip every 7 milliseconds
  for 20 repitions of the syllable 'bob'.
}
\usage{
lip
lipmarks
liptime
}
\format{
  \describe{
    \item{lip}{
      a matrix of dimension c(51, 20) giving the position of the lower
      lip every 7 milliseconds for 350 milliseconds.
    }
    \item{lipmarks}{
      a matrix of dimension c(20, 2) giving the positions of the
      'leftElbow' and 'rightElbow' in each of the 20 repetitions of the
      syllable 'bob'.
    }
    \item{liptime}{
      time in seconds from the start = seq(0, 0.35, 51) = every 7
      milliseconds.
    }
  }
}
\details{
  These are rather simple data, involving the movement of the lower lip
  while saying "bob".  There are 20 replications and 51 sampling points.
  The data are used to illustrate two techniques:  landmark registration
  and principal differental analysis.
  Principal differential analysis estimates a linear differential equation
  that can be used to describe not only the observed curves, but also a
  certain number of their derivatives.
  For a rather more elaborate example of principal differential analysis,
  see the handwriting data.

  See the \code{lip} \code{demo}.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
%\seealso{\code{\link{eval.penalty}}}
\examples{
#  See the lip demo.
}
% docclass is function
\keyword{datasets}
