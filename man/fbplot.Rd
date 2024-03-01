\name{fbplot}
\Rdversion{1.1}
\alias{fbplot}
\alias{boxplot.fd}
\alias{boxplot.fdPar}
\alias{boxplot.fdSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Functional Boxplots
}
\description{
  Produces functional boxplots or enhanced functional boxplots of the
  given functional data. It can also be used to carry out functional
  data ordering based on band depth.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fbplot(fit, x = NULL, method = "MBD", depth = NULL, plot = TRUE,
	 prob = 0.5, color = 6, outliercol = 2, barcol = 4,
	 fullout=FALSE, factor=1.5,xlim=c(1,nrow(fit)),
	 ylim=c(min(fit)-.5*diff(range(fit)),max(fit)+.5*diff(range(fit))),...)
\method{boxplot}{fd}(x, z=NULL, ...)
\method{boxplot}{fdPar}(x, z=NULL, ...)
\method{boxplot}{fdSmooth}(x, z=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
    a p-by-n functional data matrix where n is the number of curves, and
    p is defined below.
  }
  \item{x}{
    For \code{fbplot}, \code{x} is the x coordinates of curves. Defaults
    to 1:p where p is the number of x coordinates.

    For \code{boxplot.fd}, \code{boxplot.fdPar} and
    \code{boxplot.fdSmooth}, \code{x} is an object of class \code{fd},
    \code{fdPar} or \code{fdSmooth}, respectively.
  }
  \item{z}{
    The coordinate of the curves, labeled \code{x} for \code{fdplot}.
    For \code{boxplot.fd}, \code{boxplot.fdPar} and
    \code{boxplot.fdSmooth}, this cannot be \code{x}, because that would
    clash with the generic \code{boxplot(x, ...)} standard.
  }
  \item{method}{
    the method to be used to compute band depth. Can be one of "BD2",
    "MBD" or "Both" with a default of "MBD". See also details.
  }
  \item{depth}{
    a vector giving band depths of curves. If missing, band depth
    computation is conducted.
  }
  \item{plot}{
    logical. If TRUE (the default) then a functional boxplot is
    produced. If not, band depth and outliers are returned.
  }
  \item{prob}{
    a vector giving the probabilities of central regions in a decreasing
    order, then an enhanced functional boxplot is produced. Defaults to
    be 0.5 and a functional boxplot is plotted.
  }
  \item{color}{
    a vector giving the colors of central regions from light to dark for
    an enhanced functional boxplot. Defaults to be magenta for a
    functional boxplot.
  }
  \item{outliercol}{
    color of outlying curves. Defaults to be red.
  }
  \item{barcol}{
    color of bars in a functional boxplot. Defaults to be blue.
  }
  \item{fullout}{
    logical for plotting outlying curves. If FALSE (the default) then
    only the part outside the box is plotted. If TRUE, complete outlying
    curves are plotted.
  }
  \item{factor}{
    the constant factor to inflate the middle box and determine fences
    for outliers.  Defaults to be 1.5 as in a classical boxplot.
  }
  \item{xlim}{ x-axis limits }
  \item{ylim}{ y-axis limits }
  \item{...}{
    For \code{fbplot}, optional arguments for \code{plot}.

    For \code{boxplot.fd}, \code{boxplot.fdPar}, or
    \code{boxplot.fdSmooth}, optional arguments for \code{fbplot}.
  }
}
\details{
  For functional data, the band depth (BD) or modified band depth (MBD)
  allows for ordering a sample of curves from the center outwards and,
  thus, introduces a measure to define functional quantiles and the
  centrality or outlyingness of an observation. A smaller rank is
  associated with a more central position with respect to the sample
  curves.  BD usually provides many ties (curves have the same depth
  values), but MBD does not.  "BD2" uses two curves to determine a
  band. The method "Both" uses "BD2" first  and then uses "MBD" to break
  ties.  The method "Both" uses BD2 first and then uses MBD to break
  ties.  The computation is carried out by the fast algorithm proposed
  by Sun et. al. (2012).
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
  Sun, Y., Genton, M. G. and Nychka, D. (2012), "Exact fast computation
  of band depth for large functional datasets: How quickly can one
  million curves be ranked?" Stat, 1, 68-74.

  Sun, Y. and Genton, M. G. (2011), "Functional Boxplots," Journal of
  Computational and Graphical Statistics, 20, 316-334.

  Lopez-Pintado, S. and Romo, J. (2009), "On the concept of depth for
  functional data," Journal of the American Statistical Association,
  104, 718-734.

  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\author{
  Ying Sun \email{sunwards@stat.osu.edu}

  Marc G. Genton \email{marc.genton@kaust.edu.sa}
}
\examples{
##
## 1.  generate 50 random curves with some covariance structure
##     model 1 without outliers
##
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
ydata  = mu+t(L)\%*\%e

#functional boxplot
oldpar <- par(no.readonly=TRUE)
fbplot(ydata,method='MBD',ylim=c(-11,15))

# The same using boxplot.fd
boxplot.fd(ydata, method='MBD', ylim=c(-11, 15))

# same with default ylim
boxplot.fd(ydata)

##
## 2.  as an fd object
##
T      = dim(ydata)[1]
time   = seq(0,T,len=T)
ybasis = create.bspline.basis(c(0,T), 23)
Yfd    = smooth.basis(time, ydata, ybasis)$fd
boxplot(Yfd)

##
## 3.  as an fdPar object
##
Ypar <- fdPar(Yfd)
boxplot(Ypar)

##
## 4.  Smoothed version
##
Ysmooth <- smooth.fdPar(Yfd)
boxplot(Ysmooth)

##
## 5.  model 2 with outliers
##
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
par(oldpar)
}
\keyword{ ~Functional Boxplots }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
