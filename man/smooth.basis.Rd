\name{smooth.basis}
\alias{smooth.basis}
\title{Smooth Data using a Roughness Penalty }
\description{
This is the main function for smoothing data using a roughness penalty. Unlike function data2fd, which does not employ a rougness penalty, this function controls the nature and degree of smoothing by penalyzing a measure of rougness. Roughness is definable in a wide variety of ways using either derivatives or a linear differential operator. 
}
\usage{
smooth.basis(y, argvals, basisfd, wtvec=rep(1, n), Lfd=NULL, lambda=0,
             fdnames=list(NULL, dimnames(y)[2], NULL))
}
\arguments{
  \item{y }{An array containing values of curves at discrete sampling points or argument values. If the array is a matrix, the rows must correspond to argumentvalues and columns to replications, and it will be assumed that there is only one variable per observation. If Y is a three-dimensional array, the first dimension corresponds to argument values, the second to replications, and the third to variables within replications. If Y is a vector, only one replicate and variable are assumed. }
  \item{argvals }{A vector of argument values correspond to the observations in array y. }
  \item{basisfd }{A basis object. Each curve is represented by a linear combination of the basis functions defined in this object. }
  \item{wtvec }{A vector the same length as argvals containing nonnegative weight to be applied to the observations. By default these are one's. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
  \item{lambda }{A nonnegative value controlling the amount of roughness in the data. }
  \item{fdnames}{ A list of length 3 with members containing 
    a single name for the argument domain, such as 'Time' 
    a vector of names for the replications or cases 
    a name for the function, or a vector of names if there are multiple functions. 
   }
}
\details{
If lambda is zero, there is no penalty on roughness. As lambda increases, usually in logarithmic terms, the penalty on roughness increases and the fitted curves become more and more smooth. Ultimately, the curves are forced to have zero roughness in the sense of being in the null space of the differential operator Lfd. For example, a common choice of roughness penalty is the integrated square of the second derivative. This penalizes curvature. Since the second derivative of a straight line is zero, very large values of lambda will force the fit to become linear. 

It is also possible to control the amount of roughness by using a degrees of freedom measure. The value equivalent to lambda is found in the list returned by the function. On the other hand, it is possible to specify a degrees of freedom value, and then use function df2lambda to determine the equivalent value of lambda. 

One should not put complete faith in any automatic method for selecting lambda, including the GCV method. There are many reasons for this. For example, if derivatives are required, then the smoothing level that is automatically selected may give unacceptably rough derivatives. These methods are also highly sensitive to the assumption of independent errors, which is usually dubious with functional data. The best advice is to start with the value minimizing the gcv measure, and then explore lambda values a few log units up and down from this value to see what the smoothing function and its derivatives look like. The function plotfit.fd was designed for this purpose. 

An alternative to using smooth.basis is to first represent the data in a basis system with reasonably high resolution using data2fd, and then smooth the resulting functional data object using function smooth.fd. 
}
\value{
A list containing: 
fd
An object of class fd containing coefficients 
df
A degrees of freedom measure 
gcv
A measure of lack of fit discounted for df. One method for choosing the smoothing parameter lambda, called the GCV method, is to find the value which minimizes this measure. 
}
\references{
A discussion of roughness penalties can be found in Chapter 4 of Ramsay, J. O. and Silverman, B.W. (1997) Functional Data Analysis. More information can be found in recent texts on nonparametric regression. A good discussion of degrees of freedom measure can be found in Hastie, T. and Tibshirani, R. (1990) Generalized Additive Models. 
}
\note{}
\author{}
\seealso{
data2fd, plotfit.fd, smooth.fd, project.basis 
}
\examples{}
\keyword{}
