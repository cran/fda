\name{pinch}
\alias{pinch}
\non_function{}
\title{
  pinch force data
}
\description{
  151 measurements of pinch force during 20 replications  
}
\format{
  Matrix of dimension c(151, 20).  
}
\details{
  Measurements every 2 miliseconds.  
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York, p. 13, Figure 1.11,
  pp. 22-23, Figure 2.2, and p. 144, Figure 7.13.  
}
%\seealso{\code{\link{eval.penalty}}}
\examples{
plot(pinch[, 1], type="b")
}
% docclass is function
\keyword{datasets}
