\name{int2Lfd}
\alias{int2Lfd}
\title{
  Convert Integer to Linear Differential Operator
}
\description{
This function turns an integer specifying an order of a derivative into the
equivalent linear differential operator object.  It is also useful for
checking that an object is of the "Lfd" class.
}
\usage{
int2Lfd(m=0)
}
\arguments{
\item{m}{
either a nonnegative integer or a linear differential operator object.
}
}
\value{
e linear differential operator object of the "Lfd" class that is
equivalent to the integer argument.
}
% docclass is function
\keyword{smooth}
