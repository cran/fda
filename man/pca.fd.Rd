\name{pca.fd}
\alias{pca.fd}
\title{Functional Principal Components Analysis }
\description{
Functional Principal components analysis aims to display types of variation across a sample of functions. Principal components analysis is an exploratory data analysis that tends to be an early part of many projects. These modes of variation are called principal components or harmonics. This function computes these harmonics, the eigenvalues that indicate how important each mode of variation, and harmonic scores for individual functions. If the functions are multivariate, these harmonics are combined into a composite function that summarizes joint variation among the several functions that make up a multivariate functional observation. 
}
\usage{
pca.fd(fd, nharm=2, lambda=0, Lfd=2, centerfns=TRUE)
}
\arguments{
  \item{fd }{A functional data object defining a sample of curves. }
  \item{nharm }{The number of principal components or harmonics that are to be estimated for these functional data. }
  \item{lambda }{A smoothing parameter used to control the smoothness of the estimated harmonics. }
  \item{Lfd }{A linear differential operator object that defines the nature of the smoothness that is required of the harmonics. }
  \item{centerfns }{A logical variable which, if true, causes the mean function to be subtracted from each function prior to computing harmonics. That is, if true, the modes of variation are about the mean. If false, the modes of variation are about zero. }
}
\details{
}
\value{
Returns a list with the following entries: 
"harmonics" A functional data object containing the nharm harmonic, principal component, or eigenfunctions. If there is more than one variable in the fd argument, there is a harmonic corresponding to each function, and in this case the coefficient matrix has three dimensions. 
"values" The complete set of eigenvalues. The number of these is normally equal to the minimum of (1) number of basis functions defining the functional observations, or (2) the sample size. 
"scores" The principal component scores for each replication and harmonic. 
"varprop" The proportion of variance accounted for by each harmonic. 
"meanfd" meanfd The mean functional data object. 
}
\references{
Ramsay, J. O. and Silverman, B. W. (1997) Functional Data Analysis. New York: Springer. Compution of the principal component functions, eigenvalues, and principal component scores described in Chapter 6, and also incorporates the regularization concept described in Chapter 7. In Chapter 8 contains extensions to permit the analysis of joint variation in multivariate and functional data. 
}
\note{}
\author{}
\seealso{
plotPcaFd, cca.fd 
}
\examples{}
\keyword{}
