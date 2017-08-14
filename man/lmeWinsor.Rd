\name{lmeWinsor}
\alias{lmeWinsor}
\title{
  Winsorized Regression with mixed effects
}
\description{
  Clip inputs and mixed-effects predictions to (upper, lower) or to
  selected quantiles to limit wild predictions outside the training
  set.
}
\usage{
  lmeWinsor(fixed, data, random, lower=NULL, upper=NULL, trim=0,
        quantileType=7, correlation, weights, subset, method,
        na.action, control, contrasts = NULL, keep.data=TRUE,
        ...)
}
\arguments{
  \item{fixed}{
    a two-sided linear formula object describing the fixed-effects
    part of the model, with the response on the left of a '~'
    operator and the terms, separated by '+' operators, on the right.
    The left hand side of 'formula' must be a single vector in 
    'data', untransformed.
  }
  \item{data}{
    an optional data frame containing the variables named in 'fixed',
    'random', 'correlation', 'weights', and 'subset'.  By default the
    variables are taken from the environment from which 
    \link[nlme]{lme} is  called.
  }
  \item{random}{
    a random- / mixed-effects specification, as described with
    \link[nlme]{lme}.

    NOTE:  Unlike \link[nlme]{lme}, 'random' must be provided;  it 
    can not be inferred from 'data'.
  }
  \item{lower, upper}{
    optional numeric vectors with names matching columns of 'data'
    giving limits on the ranges of predictors and predictions:  If
    present, values below 'lower' will be increased to 'lower', and
    values above 'upper' will be decreased to 'upper'.  If absent, 
    these limit(s) will be inferred from quantile(..., prob=c(trim, 
    1-trim), na.rm=TRUE, type=quantileType).
  }
  \item{trim}{
    the fraction (0 to 0.5) of observations to be considered outside 
    the range of the data in determining limits not specified in
    'lower' and 'upper'.

    NOTES:

    (1) trim>0 with a singular fit may give an error.  In such cases,
    fix the singularity and retry.

    (2) trim = 0.5 should NOT be used except to check the algorithm,
    because it trims everything to the median, thereby providing zero
    leverage for estimating a regression.

    (3) The current algorithm does does NOT adjust any of the 
    variance parameter estimates to account for predictions outside
    'lower' and 'upper'.  This will have no effect for trim = 0 or
    trim otherwise so small that there are not predictions outside
    'lower' and 'upper'.  However, for more substantive trimming,
    this could be an issue. This is different from \link{lmWinsor}.
  }
  \item{quantileType}{
    an integer between 1 and 9 selecting one of the nine quantile
    algorithms to be used with 'trim' to determine limits not
    provided with 'lower' and 'upper'.
  }
  \item{correlation}{
    an optional correlation structure, as described with
    \link[nlme]{lme}.
  }
  \item{weights}{
    an optional heteroscedasticity structure, as described with
    \link[nlme]{lme}.
  }
  \item{ subset }{
    an optional vector specifying a subset of observations to be 
    used in the fitting process, as described with \link[nlme]{lme}.
  }
  \item{method}{
    a character string.  If '"REML"' the model is fit by maximizing 
    the restricted log-likelihood.  If '"ML"' the log-likelihood is
    maximized.  Defaults to '"REML"'.
  }
  \item{ na.action }{
    a function that indicates what should happen when the data 
    contain 'NA's.  The default action ('na.fail') causes 'lme' to 
    print an error message and terminate if there are any incomplete
    observations.
  }
  \item{control}{
    a list of control values for the estimation algorithm to replace 
    the default values returned by the function
    \link[nlme]{lmeControl}. Defaults to an empty list.

    NOTE:  Other control parameters such as 'singular.ok' as
    documented in \link[nlme]{glsControl} may also work, but should
    be used with caution.
  }
  \item{ contrasts }{
    an optional list. See the 'contrasts.arg' of
    'model.matrix.default'.
  }
  \item{keep.data}{
    logical: should the 'data' argument (if supplied and a 
    data.frame) be saved as part of the model object?
  }
  \item{\dots}{
    additional arguments to be passed to the low level regression
    fitting functions;  see \link[nlme]{lme}.
  }
}
\details{
  1.  Identify inputs and outputs as follows:

  1.1.  mdly <- mdlx <- fixed; mdly[[3]] <- NULL; mdlx[[2]] <- NULL;

  1.2.  xNames <- c(all.vars(mdlx), all.vars(random)).

  1.3.  yNames <- all.vars(mdly).  Give an error if
  as.character(mdly[[2]]) != yNames.

  2.  Do 'lower' and 'upper' contain limits for all numeric columns 
  of 'data'?  Create limits to fill any missing.

  3.  clipData = data with all xNames clipped to (lower, upper).

  4.  fit0 <- lme(...)

  5.  Add components lower and upper to fit0 and convert it to class
  c('lmeWinsor', 'lme').

  6.  Clip any stored predictions at the Winsor limits for 'y'.

  NOTE:  This is different from \link{lmWinsor}, which uses quadratic
  programming with predictions outside limits, transferring extreme
  points one at a time to constraints that force the unWinsorized
  predictions for those points to be at least as extreme as the
  limits.
}
\value{
  an object of class c('lmeWinsor', 'lme') with 'lower', 'upper', and
  'message' components in addition to the standard 'lm' components.  
  The 'message' is a list with its first component being either 'all
  predictions inside limits' or 'predictions outside limits'.  In the
  latter case, there rest of the list summarizes how many and which
  points have predictions outside limits.
}
\author{ Spencer Graves }
\seealso{
  \code{\link{lmWinsor}}
  \code{\link{predict.lmeWinsor}}
  \code{\link[nlme]{lme}}
  \code{\link{quantile}}
}
\examples{
requireNamespace('nlme')
fm1w <- lmeWinsor(distance ~ age, data = nlme::Orthodont,
                 random=~age|Subject)
fm1w.1 <- lmeWinsor(distance ~ age, data = nlme::Orthodont,
                 random=~age|Subject, trim=0.1)
}
\keyword{ models }
