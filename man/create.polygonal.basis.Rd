\name{create.polygonal.basis}
\alias{create.polygonal.basis}
\title{Create a Polygonal Basis }
\description{
A basis is set up for constructing polygonal lines, consisting of straight line segments that join together. 
}
\usage{
create.polygonal.basis(argvals=NULL)
}
\arguments{
  \item{argvals }{A strictly increasing vector of argument values at which line segments join to form a polygonal line. }
}
\details{
The actual basis functions consist of triangles, each with its apex over an argument value. Note that in effect the polygonal basis is identical to a B-spline basis of order 2 and a knot or break value at each argument value. The range of the polygonal basis is set to the interval defined by the smallest and largest argument values. 
}
\value{
A basis object with the type "polyg". 
}
\references{}
\note{}
\author{}
\seealso{
basis.fd. create.basis, create.bspline.basis 
}
\examples{}
\keyword{}
