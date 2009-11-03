\name{fbplot}
\Rdversion{1.1}
\alias{fbplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional Boxplots
}
\description{
Produces functional boxplots or enhanced functional boxplots of the given functional data.
It can also be used to carry out functional data ordering based on band depth. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fbplot(fit, x = NULL, method = "MBD", depth = NULL, plot = TRUE, 
	 prob = 0.5, color = 6, outliercol = 2, barcol = 4, fullout=FALSE, 
   factor = 1.5, ...)	 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
  either a functional data object or 
  a p-by-n data matrix where n is the number of curves, and p is the number 
  of evaluation points. 
%%     ~~Describe \code{fit} here~~
}
  \item{x}{
  evaluation points of the curves. Defaults to 101 equally spaced
  points on the range of \code{fit} if \code{fit} is a functional data object
  or to 1:p where p is the number of x coordinates if \code{fit} is a data matrix.
%%     ~~Describe \code{x} here~~
}
  \item{method}{
  the method to be used to compute band depth. Can be one of "BD2", "BD3", "MBD" or "Both" with a default of "MBD". See also details.
%%     ~~Describe \code{method} here~~
}
  \item{depth}{
  a vector giving band depths of curves. If missing, band depth computation is conducted.
%%     ~~Describe \code{depth} here~~
}
  \item{plot}{
  logical. If TRUE (the default) then a functional boxplot is produced. If not, band depth and outliers are returned.
%%     ~~Describe \code{plot} here~~
}
  \item{prob}{
  a vector giving the probabilities of central regions in a decreasing order, then an enhanced functional boxplot is produced. Defaults to be 0.5 and a functional boxplot is plotted.
%%     ~~Describe \code{prob} here~~
}
  \item{color}{
  a vector giving the colors of central regions from light to dark for an enhanced functional boxplot. Defaults to be magenta for a functional boxplot.
%%     ~~Describe \code{color} here~~
}
  \item{outliercol}{
  color of outlying curves. Defaults to be red.
%%     ~~Describe \code{outliercol} here~~
}
  \item{barcol}{
  color of bars in a functional boxplot. Defaults to be blue.
%%     ~~Describe \code{barcol} here~~
}
 \item{fullout}{
  logical for plotting outlying curves. If FALSE (the default) then only the part  outside the box is plotted. If TRUE, complete outling curves are plotted.
%%     ~~Describe \code{fullout} here~~
}
\item{factor}{
  the constant factor to inflate the middle box and determine fences for outliers. Defaults to be 1.5 as in a classical boxplot.
}
%%     ~~Describe \code{factor} here~~
\item{...}{ additional plotting arguments}
}
\details{
For functional data, the band depth (BD) or modifed band depth (MBD) allows for ordering a sample of curves from the center outwards 
and, thus, introduces a measure to define functional quantiles and the centrality or outlyingness of an
observation. A smaller rank is associated with a more central position with respect to the sample curves. "BD2" uses two curves to determine a band and "BD3" uses three curves.
BD usually provides many ties (curves have the same depth values), but MBD does not. The method "Both" uses BD2 first and then uses MBD to break ties.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{depth }{band depths of given curves.}
\item{outpoint }{column indices of detected outliers.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Sun, Y., and Genton, M. G. (2011), "Functional boxplots,"
<http://www.stat.tamu.edu/%7Egenton/2011.SG.JCGS.pdf> /Journal of
Computational and Graphical Statistics/, to appear.
\cr
\cr Lopez-Pintado, S. and Romo, J. (2009), "On the concept of depth for functional data," Journal
of the American Statistical Association, 104, 718-734.
%% ~put references to the literature/web site here ~
}
\author{
Ying Sun \email{sunwards@stat.tamu.edu}
\cr
\cr Marc G. Genton \email{genton@stat.tamu.edu}
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#generate 50 random curves with some covariance structure
#model 1 without outliers
cov.fun=function(d,k,c,mu){
        k*exp(-c*d^mu)
}
n=50
p=30
t=seq(0,1,len=p)
d=dist(t,upper=TRUE,diag=TRUE)
d.matrix=as.matrix(d)
#covariance function in time
t.cov=cov.fun(d.matrix,1,1,1)
# Cholesky Decomposition
L=chol(t.cov)
mu=4*t
e=matrix(rnorm(n*p),p,n)
y=mu+t(L)%*%e

#functional boxplot
fbplot(y,method='MBD',ylim=c(-11,15))

#model 2 with outliers
#magnitude
k=6
#randomly introduce outliers
C=rbinom(n,1,0.1)
s=2*rbinom(n,1,0.5)-1
cs.m=matrix(C*s,p,n,byrow=TRUE)

e=matrix(rnorm(n*p),p,n)
y=mu+t(L)\%*\%e+k*cs.m

#functional boxplot
fbplot(y,method='MBD',ylim=c(-11,15))

}
\keyword{ ~Functional Boxplots }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
