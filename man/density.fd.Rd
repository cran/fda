\name{density.fd}
\alias{density.fd}
\title{
  Compute a Probability Density Function
}
\description{
  Like the regular S-PLUS function \code{density}, this function
  computes a probability density function for a sample of values of a
  random variable.  However, in this case the density function is
  defined by a functional parameter object \code{WfdParobj} along with a
  normalizing constant \code{C}.

  The density function $p(x)$ has the form \code{p(x) = C exp[W(x)]}
  where function $W(x)$ is defined by the functional data object
  \code{WfdParobj}.
}
\usage{
\method{density}{fd}(x, WfdParobj, conv=0.0001, iterlim=20,
           active=1:nbasis, dbglev=0, ...)
}
\arguments{
\item{x}{
a set observations, which may be one of two forms:
\enumerate{
	\item a vector of observations $x_i$
	\item a two-column matrix, with the observations $x_i$ in the
	first column, and frequencies $f_i$ in the second.
}
The first option corresponds to all $f_i = 1$.
}
\item{WfdParobj}{
a functional parameter object specifying the initial
value, basis object, roughness penalty and smoothing
parameter defining function $W(t).$
}
\item{conv}{
a positive constant defining the convergence criterion.
}
\item{iterlim}{
the maximum number of iterations allowed.
}
\item{active}{
a logical vector of length equal to the number of coefficients
defining \code{Wfdobj}. If an entry is TRUE, the corresponding
coefficient is estimated, and if FALSE, it is held at the value defining the
argument \code{Wfdobj}.  Normally the first coefficient is set to 0
and not estimated, since it is assumed that $W(0) = 0$.
}
\item{dbglev}{
either 0, 1, or 2.  This controls the amount information printed out on
each iteration, with 0 implying no output, 1 intermediate output level,
and 2 full output.  If levels 1 and 2 are used, it is helpful to
turn off the output buffering option in S-PLUS.
}
\item{...}{Other arguments to match the generic function 'density'}
}
\value{
  a named list of length 4 containing:
  \item{Wfdobj}{
    a functional data object defining function $W(x)$ that that
    optimizes the fit to the data of the monotone function that it defines.
  }
  \item{C}{
    the normalizing constant.
  }
  \item{Flist}{
    a named list containing three results for the final converged solution:
    (1)
    \bold{f}: the optimal function value being minimized,
    (2)
    \bold{grad}: the gradient vector at the optimal solution,   and
    (3)
    \bold{norm}: the norm of the gradient vector at the optimal solution.
  }
  \item{iternum}{
    the number of iterations.
  }
  \item{iterhist}{
    a \code{iternum+1} by 5 matrix containing the iteration
    history.
  }
}
\details{
  The goal of the function is provide a smooth density function
  estimate that approaches some target density by an amount that is
  controlled by the linear differential operator \code{Lfdobj} and
  the penalty parameter. For example, if the second derivative of
  $W(t)$ is penalized heavily, this will force the function to
  approach a straight line, which in turn will force the density function
  itself to be nearly normal or Gaussian.  Similarly, to each textbook
  density function there corresponds a $W(t)$, and to each of these
  in turn their corresponds a linear differential operator that will, when
  apply to $W(t)$, produce zero as a result.
  To plot the density function or to evaluate it, evaluate \code{Wfdobj},
  exponentiate the resulting vector, and then divide by the normalizing
  constant \code{C}.
}
\seealso{
  \code{\link{intensity.fd}}
}
\examples{

#  set up range for density
rangeval <- c(-3,3)
#  set up some standard normal data
x <- rnorm(50)
#  make sure values within the range
x[x < -3] <- -2.99
x[x >  3] <-  2.99
#  set up basis for W(x)
basisobj <- create.bspline.basis(rangeval, 11)
#  set up initial value for Wfdobj
Wfd0 <- fd(matrix(0,11,1), basisobj)
WfdParobj <- fdPar(Wfd0)
#  estimate density
denslist <- density.fd(x, WfdParobj)
#  plot density
oldpar <- par(no.readonly=TRUE)
xval <- seq(-3,3,.2)
wval <- eval.fd(xval, denslist$Wfdobj)
pval <- exp(wval)/denslist$C
plot(xval, pval, type="l", ylim=c(0,0.4))
points(x,rep(0,50))
par(oldpar)
}
% docclass is function
\keyword{smooth}
