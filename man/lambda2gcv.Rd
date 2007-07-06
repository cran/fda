\name{lambda2gcv}
\alias{lambda2gcv}
\title{
  Compute GCV Criterion
}
\description{
The generalized cross-validation or GCV criterion is often
used to select an appropriate smoothing parameter value, by
finding the smoothing parameter that minimizes GCV.  This
function locates that value.
}
\usage{
%lambda2gcv(log10lambda, argvals, y, fdParobj,
%           wtvec=rep(1, length(y)), dffactor=1)
lambda2gcv(log10lambda, argvals, y, fdParobj,
           wtvec=rep(1, length(y)))
}
\arguments{
\item{log10lambda}{
the logarithm (base 10) of the smoothing parameter
}
\item{argvals}{
a vector of argument values.
}
\item{y}{
the data to be smoothed.
}
\item{fdParobj}{
a functional parameter object defining the smooth.
}
\item{wtvec}{
a weight vector used in the smoothing.
}
%\item{dffactor}{a factor used to define the GCV criterion.  See
%\code{smooth.basis} for details.}
}
\value{
the value of the GCV criterion.
}
% docclass is function
\keyword{smooth}
