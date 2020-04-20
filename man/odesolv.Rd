\name{odesolv}
\alias{odesolv}
\title{
  Numerical Solution mth Order Differential Equation System
}
\description{
  The system of differential equations is linear, with possibly
  time-varying coefficient functions.  The numerical solution is
  computed with the Runge-Kutta method.
}
\usage{
odesolv(bwtlist, ystart=diag(rep(1,norder)),
        h0=width/100, hmin=width*1e-10, hmax=width*0.5,
        EPS=1e-4, MAXSTP=1000)
}
\arguments{
\item{bwtlist}{
  a list whose members are functional parameter objects defining the
  weight functions for the linear differential equation.
}
\item{ystart}{
  a vector of initial values for the equations.  These are the values at
  time 0 of the solution and its first m - 1 derivatives.
}
\item{h0}{
a positive initial step size.
}
\item{hmin}{
the minimum allowable step size.
}
\item{hmax}{
the maximum allowable step size.
}
\item{EPS}{
a convergence criterion.
}
\item{MAXSTP}{
the maximum number of steps allowed.
}
}
\value{
a named list of length 2 containing

\item{tp}{
a vector of time values at which the system is evaluated
}
\item{yp}{
a matrix of variable values corresponding to \code{tp}.
}
}
\details{
This function is required to compute a set of solutions of an
estimated linear differential equation in order compute a fit
to the data that solves the equation.  Such a fit will be a
linear combinations of m independent solutions.
}
\seealso{
  \code{\link{pda.fd}}.  For new applications, users are encouraged to
  consider \code{\link[deSolve]{deSolve}}.  The \code{deSolve} package
  provides general solvers for ordinary and partial differential
  equations, as well as differential algebraic equations and delay
  differential equations.
}
\examples{
#See the analyses of the lip data.
}
% docclass is function
\keyword{smooth}
