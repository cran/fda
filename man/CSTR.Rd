\name{CSTR}
\alias{CSTR}
\alias{CSTR2in}
\alias{CSTR2}
\alias{CSTRfn}
\alias{CSTRfitLS}
\alias{CSTRres}
\alias{CSTRsse}
\title{
  Continuously Stirred Tank Reactor
}
\description{
  Functions for solving the Continuously Stirred Tank Reactor
  (CSTR) Ordinary Differential Equations (ODEs).  A solution for
  observations where metrology error is assumed to be negligible can be
  obtained via lsoda(y, Time, CSTR2, parms);  CSTR2 calls CSTR2in.  When
  metrology error can not be ignored, use CSTRfn (which calls
  CSTRfitLS).  To estimate parameters in the CSTR differential equation
  system (kref, EoverR, a, and / or b), pass CSTRres to nls.  If nls
  fails to converge, first use optim or nlminb with CSTRsse, then pass
  the estimates to nls.
}
\usage{
CSTR2in(Time, condition =
   c('all.cool.step', 'all.hot.step', 'all.hot.ramp', 'all.cool.ramp',
     'Tc.hot.exponential', 'Tc.cool.exponential', 'Tc.hot.ramp',
     'Tc.cool.ramp', 'Tc.hot.step', 'Tc.cool.step'),
   tau=1)
CSTR2(Time, y, parms)

CSTRfitLS(coef, datstruct, fitstruct, lambda, gradwrd=FALSE)
CSTRfn(parvec, datstruct, fitstruct, CSTRbasis, lambda, gradwrd=TRUE)
CSTRres(kref=NULL, EoverR=NULL, a=NULL, b=NULL,
        datstruct, fitstruct, CSTRbasis, lambda, gradwrd=FALSE)
CSTRsse(par, datstruct, fitstruct, CSTRbasis, lambda)

}
\arguments{
  \item{Time}{
    The time(s) for which computation(s) are desired
  }
  \item{condition}{
    a character string with the name of one of ten preprogrammed input
    scenarios.
  }
  \item{ tau }{
    time for exponential decay of exp(-1) under condition =
    'Tc.hot.exponential' or 'Tc.cool.exponential';  ignored for other
    values of 'condition'.
  }
  \item{y}{
    Either a vector of length 2 or a matrix with 2 columns giving the
    observation(s) on Concentration and Temperature for which
    computation(s) are desired
  }

  \item{parms}{
    a list of CSTR model parameters passed via the lsoda 'parms'
    argument.  This list consists of the following 3 components:

    \itemize{
      \item{fitstruct}{
	a list with 12 components describing the structure for fitting.
	This is the same as the 'fitstruct' argument of 'CSTRfitLS' and
	'CSTRfn' without the 'fit' component;  see below.
      }

      \item{condition}{
	a character string identifying the inputs to the simulation.
	Currently, any of the following are accepted:  'all.cool.step',
	'all.hot.step', 'all.hot.ramp', 'all.cool.ramp',
	'Tc.hot.exponential', 'Tc.cool.exponential', 'Tc.hot.ramp',
	'Tc.cool.ramp', 'Tc.hot.step', or 'Tc.cool.step'.
      }
      \item{Tlim}{
	end time for the computations.
      }
    }
  }
  \item{coef}{
    a matrix with one row for each basis function in fitstruct and
    columns c("Conc", "Temp") or a vector form of such a matrix.
  }
  \item{datstruct}{
    a list describing the structure of the data.  CSTRfitLS uses the
    following components:

    \itemize{
      \item{basismat, Dbasismat}{
	basis coefficent matrices with one row for each observation and
	one column for each basis vector.  These are typically produced
	by code something like the following:

	basismat <- eval.basis(Time, CSTRbasis)

	Dbasismat <- eval.basis(Time, CSTRbasis, 1)
      }
      \item{Cwt, Twt}{
	scalar variances of 'fd' functional data objects for
	Concentration and Temperature used to place the two series on
	comparable scales.
      }

      \item{y}{
	a matrix with 2 columns for the observed 'Conc' and 'Temp'.
      }
      \item{quadbasismat, Dquadbasismat}{
	basis coefficient matrices with one row for each quadrature
	point and one column for each basis vector.  These are typically
	produced by code something like the following:

	quadbasismat <- eval.basis(quadpts, CSTRbasis)

	Dquadbasismat <- eval.basis(quadpts, CSTRbasis, 1)
      }

      \item{Fc, F., CA0, T0, Tc}{
	input series for CSTRfitLS and CSTRfn as the output list
	produced by CSTR2in.
      }

      \item{quadpts}{
	Quadrature points created by 'quadset' and stored in
	CSTRbasis[["quadvals"]][, "quadpts"].
      }
      \item{quadwts}{
	Quadrature weights created by 'quadset' and stored in
	CSTRbasis[["quadvals"]][, "quadpts"].
      }
    }
  }
  \item{fitstruct}{
    a list with 14 components:

    \itemize{
      \item{V}{
	volume in cubic meters
      }
      \item{Cp}{
	concentration in cal/(g.K) for computing betaTC and betaTT;  see
	details below.
      }
      \item{rho}{
	density in grams per cubic meter
      }
      \item{delH}{
	cal/kmol
      }
      \item{Cpc}{
	concentration in cal/(g.K) used for computing alpha;  see
	details below.
      }
      \item{Tref}{
	reference temperature.
      }
      \item{kref}{
	reference value
      }
      \item{EoverR}{
	E/R in units of K/1e4
      }
      \item{a}{
	scale factor for Fco in alpha;  see details below.
      }
      \item{b}{
	power of Fco in alpha;  see details below.
      }
      \item{Tcin}{
	Tc input temperature vector.
      }
      \item{fit}{
	logical vector of length 2 indicating whether Contentration or
	Temperature or both are considered to be observed and used for
	parameter estimation.
      }
      \item{coef0}{
	data.frame(Conc = Cfdsmth[["coef"]], Temp = Tfdsmth[["coef"]]),
	where Cfdsmth and Tfdsmth are the objects returned by
	smooth.basis when applied to the observations on Conc and Temp,
	respectively.
      }
      \item{estimate}{
	logical vector of length 4 indicating which of kref, EoverR, a
	and b are taken from 'parvec';  all others are taken from
	'fitstruct'.
      }
    }
  }
  \item{lambda}{
    a 2-vector of rate parameters 'lambdaC' and 'lambdaT'.
  }
  \item{gradwrd}{
    a logical scalar TRUE if the gradient is to be returned as well as
    the residuals matrix.
  }
  \item{parvec, par}{
    initial values for the parameters specified by fitstruct[[
    "estimate"]] to be estimated.
  }
  \item{CSTRbasis}{
    Quadrature basis returned by 'quadset'.
  }
  \item{kref, EoverR, a, b}{
    the kref, EoverR, a, and b coefficients of the CSTR model as
    individual arguments of CSTRres to support using 'nls' with the CSTR
    model.  Those actually provided by name will be estimated;  the
    others will be taken from '.fitstruct';  see details.
  }
}
\details{
  Ramsay et al. (2007) considers the following differential equation
  system for a continuously stirred tank reactor (CSTR):

  dC/dt = (-betaCC(T, F.in)*C + F.in*C.in)

  dT/dt = (-betaTT(Fcvec, F.in)*T + betaTC(T, F.in)*C +
  alpha(Fcvec)*T.co)

  where

  betaCC(T, F.in) = kref*exp(-1e4*EoverR*(1/T - 1/Tref)) + F.in

  betaTT(Fcvec, F.in) = alpha(Fcvec) + F.in

  betaTC(T, F.in) = (-delH/(rho*Cp))*betaCC(T, F.in)

  \deqn{
    alpha(Fcvec) = (a*Fcvec^(b+1) / (K1*(Fcvec + K2*Fcvec^b)))
  }

  K1 = V*rho*Cp

  K2 = 1/(2*rhoc*Cpc)

  The four functions CSTR2in, CSTR2, CSTRfitLS, and CSTRfn compute
  coefficients of basis vectors for two different solutions to this set
  of differential equations.  Functions CSTR2in and CSTR2 work with
  'lsoda' to provide a solution to this system of equations.  Functions
  CSTSRitLS and CSTRfn are used to estimate parameters to fit this
  differential equation system to noisy data.  These solutions are
  conditioned on specified values for kref, EoverR, a, and b.  The other
  function, CSTRres, support estimation of these parameters using
  'nls'.

  CSTR2in translates a character string 'condition' into a data.frame
  containing system inputs for which the reaction of the system is
  desired.  CSTR2 calls CSTR2in and then computes the corresponding
  predicted first derivatives of CSTR system outputs according to the
  right hand side of the system equations.  CSTR2 can be called by
  'lsoda' in the 'deSolve' package to actually solve the system of
  equations.  To solve the CSTR equations for another set of inputs, the
  easiest modification might be to change CSTR2in to return the desired
  inputs.  Another alternative would be to add an argument
  'input.data.frame' that would be used in place of CSTR2in when
  present.

  CSTRfitLS computes standardized residuals for systems outputs Conc,
  Temp or both as specified by fitstruct[["fit"]], a logical vector of
  length 2.  The standardization is sqrt(datstruct[["Cwt"]]) and / or
  sqrt(datstruct[["Twt"]]) for Conc and Temp, respectively.  CSTRfitLS
  also returns standardized deviations from the predicted first
  derivatives for Conc and Temp.

  CSTRfn uses a Gauss-Newton optimization to estimates the coefficients
  of CSTRbasis to minimize the weighted sum of squares of residuals
  returned by CSTRfitLS.

  CSTRres provides an interface between 'nls' and 'CSTRfn'.  It gets the
  parameters to be estimated via the official function arguments, kref,
  EoverR, a, and / or b.  The subset of these paramters to estimate must
  be specified both directly in the function call to 'nls' and
  indirectly via fitstruct[["estimate"]].  CSTRres gets the other CSTRfn
  arguments (datstruct, fitstruct, CSTRbasis, and lambda) via the 'data'
  argument of 'nls'.

  CSTRsse computes sum of squares of residuals for use with optim or
  nlminb.
}
\value{
  \item{CSTR2in}{
    a matrix with number of rows = length(Time) and columns for F., CA0,
    T0, Tcin, and Fc.  This gives the inputs to the CSTR simulation for
    the chosen 'condition'.
  }

  \item{CSTR2}{
    a list with one component being a matrix with number of rows =
    length(tobs) and 2 columns giving the first derivatives of Conc and
    Temp according to the right hand side of the differential equation.
    CSTR2 calls CSTR2in to get its inputs.
  }

  \item{CSTRfitLS}{
    a list with one or two components as follows:

    \itemize{
      \item{res}{
	a list with two components

	Sres = a matrix giving the residuals between observed and
	predicted datstruct[["y"]] divided by sqrt(datstruct[[c("Cwt",
	"Twt")]]) so the result is dimensionless.  dim(Sres) =
	dim(datstruct[["y"]]).  Thus, if datstruct[["y"]] has only one
	column, 'Sres' has only one column.

	Lres = a matrix with two columns giving the difference between
	left and right hand sides of the CSTR differential equation at
	all the quadrature points.  dim(Lres) = c(nquad, 2).
      }
      \item{Dres}{
	If gradwrd=TRUE, a list with two components:

	DSres = a matrix with one row for each element of res[["Sres"]]
	and two columns for each basis function.

	DLres = a matrix with two rows for each quadrature point and two
	columns for each basis function.

	If gradwrd=FALSE, this component is not present.
      }
    }
  }

  \item{CSTRfn}{
    a list with five components:

    \itemize{
      \item{res}{
	the 'res' component of the final 'CSTRfitLS' object reformatted
	with its component Sres first followed by Lres, using
	with(CSTRfitLS(...)[["res"]], c(Sres, Lres)).
      }
      \item{Dres}{
	one of two very different gradient matrices depending on the
	value of 'gradwrd'.

	If gradwrd = TRUE, Dres is a matrix with one row for each
	observation value to match and one column for each parameter
	taken from 'parvec' per fitstruct[["estimate"]].  Also, if
	fitstruct[["fit"]] = c(1,1), CSTRfn tries to  match both
	Concentration and Temperature, and rows corresponding to
	Concentration come first following by rows corresponding to
	Temperature.

	If gradwrd = FALSE, this is the 'Dres' component of the final
	'CSTRfitLS' object reformatted as follows:

        Dres <- with(CSTRfitLS(...)[["Dres"]], rbind(DSres, DLres))
      }
      \item{fitstruct}{
	a list components matching the 'fitstruct' input, with
	coefficients estimated replaced by their initial values from
	parvec and with coef0 replace by its final estimate.
      }
      \item{df}{
	estimated degrees of freedom as the trace of the appropriate
	matrix.
      }
      \item{gcv}{
	the Generalized cross validation estimate of the mean square
	error, as discussed in Ramsay and Silverman (2006, sec. 5.4).
      }
    }
  }
  \item{CSTRres}{
    the 'res' component of CSTRfd(...) as a column vector.  This allows
    us to use 'nls' with the CSTR model.  This can be especially useful
    as 'nls' has several helper functions to facilitate evaluating
    goodness of fit and and uncertainty in parameter estimates.
  }
  \item{CSTRsse}{
    sum(res*res) from CSTRfd(...).  This allows us to use 'optim' or
    'nlminb' with the CSTR model.  This can also be used to obtain
    starting values for 'nls' in cases where 'nls' fails to converge
    from the initiall provided starting values.  Apart from 'par', the
    other arguments 'datstruct', 'fitstruct', 'CSTRbasis', and 'lambda',
    must be passed via '...' in 'optim' or 'nlminb'.
  }
}
\references{
  Ramsay, J. O., Hooker, G., Cao, J. and Campbell, D. (2007) Parameter
  estimation for differential equations: A generalized smoothing
  approach (with discussion). \emph{Journal of the Royal Statistical
    Society}, Series B, 69, 741-796.

  Ramsay, J. O., and Silverman, B. W. (2006) \emph{Functional Data
    Analysis}, 2nd ed., Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}

\seealso{
  \code{\link[deSolve]{lsoda}}
  \code{\link{nls}}
}
\examples{
###
###
### 1.  lsoda(y, times, func=CSTR2, parms=...)
###
###
#  The system of two nonlinear equations has five forcing or
#  input functions.
#  These equations are taken from
#  Marlin, T. E. (2000) Process Control, 2nd Edition, McGraw Hill,
#  pages 899-902.
##
##  Set up the problem
##
fitstruct <- list(V    = 1.0,#  volume in cubic meters
                  Cp   = 1.0,#  concentration in cal/(g.K)
                  rho  = 1.0,#  density in grams per cubic meter
                  delH = -130.0,# cal/kmol
                  Cpc  = 1.0,#  concentration in cal/(g.K)
                  rhoc = 1.0,#  cal/kmol
                  Tref = 350)#  reference temperature
#  store true values of known parameters
EoverRtru = 0.83301#   E/R in units K/1e4
kreftru   = 0.4610 #   reference value
atru      = 1.678#     a in units (cal/min)/K/1e6
btru      = 0.5#       dimensionless exponent

#% enter these parameter values into fitstruct

fitstruct[["kref"]]   = kreftru#
fitstruct[["EoverR"]] = EoverRtru#  kref = 0.4610
fitstruct[["a"]]      = atru#       a in units (cal/min)/K/1e6
fitstruct[["b"]]      = btru#       dimensionless exponent

Tlim  = 64#    reaction observed over interval [0, Tlim]
delta = 1/12#  observe every five seconds
tspan = seq(0, Tlim, delta)#

coolStepInput <- CSTR2in(tspan, 'all.cool.step')

#  set constants for ODE solver

#  cool condition solution
#  initial conditions

Cinit.cool = 1.5965#  initial concentration in kmol per cubic meter
Tinit.cool = 341.3754# initial temperature in deg K
yinit = c(Conc = Cinit.cool, Temp=Tinit.cool)

#  load cool input into fitstruct

fitstruct[["Tcin"]] = coolStepInput[, "Tcin"];

#  solve  differential equation with true parameter values

if (require(deSolve)) {
coolStepSoln <- lsoda(y=yinit, times=tspan, func=CSTR2,
  parms=list(fitstruct=fitstruct, condition='all.cool.step', Tlim=Tlim) )
}
###
###
### 2.  CSTRfn
###
###

# See the script in '~R\library\fda\scripts\CSTR\CSTR_demo.R'
#  for more examples.

}
% docclass is function
\keyword{smooth}

