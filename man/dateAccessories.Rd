\name{monthAccessories}
\alias{monthAccessories}
\alias{dayOfYear}
\alias{day.5}
\alias{daysPerMonth}
\alias{monthEnd}
\alias{monthEnd.5}
\alias{monthBegin.5}
\alias{monthMid}
\alias{monthLetters}
\alias{weeks}
\title{
  Numeric and character vectors to facilitate working with dates 
}
\description{
  Numeric and character vectors to simplify functional data computations
  and plotting involving dates.  
}
\format{

  \item{dayOfYear}{
    a numeric vector = 1:365 
  }
  \item{day.5 }{
    a numeric vector = dayOfYear-0.5 = 0.5, 1.5, ..., 364.5 
  }
  \item{daysPerMonth}{
    a numeric vector of the days in each month (ignoring leap years)
    with names = month.abb
  }
  \item{monthEnd}{
    a numeric vector of cumsum(daysPerMonth) with names = month.abb
  }
  \item{monthEnd.5}{
    a numeric vector of the middle of the last day of each month with
    names = month.abb = c(Jan=30.5, Feb=58.5, ..., Dec=364.5)
  }
  \item{monthBegin.5}{
    a numeric vector of the middle of the first day of each month with
    names - month.abb = c(Jan=0.5, Feb=31.5, ..., Dec=334.5) 
  }
  \item{monthMid}{
    a numeric vector of the middle of the month = (monthBegin.5 +
    monthEnd.5)/2  
  }
  \item{monthLetters}{
    A character vector of c("j", "FALSE", "m", "A", "M", "J", "J", "A", "S",
    "O", "N", "D"), with 'month.abb' as the names.  
  }
  \item{weeks}{
    a numeric vector of length 53 marking 52 periods of approximately 7
    days each throughout the year = c(0, 365/52, ..., 365)
  }
}
\details{
  Miscellaneous vectors often used in 'fda' scripts.  
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York, pp. 5, 47-53.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York

}
\seealso{
  \code{\link{month.abb}} 
}
\examples{



}
% docclass is function
\keyword{datasets}

