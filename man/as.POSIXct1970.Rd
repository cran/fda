\name{as.POSIXct1970}
\alias{as.POSIXct1970}
\title{
  \code{as.POXIXct} for number of seconds since the start of 1970.
}
\description{
  \code{as.POSIXct.numeric} requires origin to be specified.  This assumes that
  it is the start of 1970.
}
\usage{as.POSIXct1970(x, tz="GMT", ...)}
\arguments{
  \item{x}{
    a numeric vector of times in seconds since the start of 1970.  (If
    \code{x} is not numeric, call \code{as.POSIXct}.)
  }
  \item{tz}{
    A timezone specification to be used for the conversion, if one is
    required. System-specific (see time zones), but "" is the current
    timezone, and "GMT" is UTC (Universal Time, Coordinated).
  }
  \item{\dots}{optional arguments to pass to as.POSIXct.}
}
\details{
  o1970 <- strptime('1970-01-01', '%Y-%m-%d', tz='GMT'),
  o1970. <- as.POSIXct(o1970),
  as.POSIXct(x, origin=o1970.)
}
\value{
  Returns a vector of class \code{POSIXct}.
}
\author{
  Spencer Graves
}
\seealso{
  \code{\link{as.POSIXct}},
  \code{\link{ISOdate}},
  \code{\link{strptime}}
}
\examples{
sec <- c(0, 1, 60, 3600, 24*3600, 31*24*3600, 365*24*3600)
Sec <- as.POSIXct1970(sec)
\dontshow{stopifnot(}
all.equal(sec, as.numeric(Sec))
\dontshow{)}
}
\keyword{manip}
