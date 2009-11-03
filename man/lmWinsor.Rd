\name{lmWinsor}
\alias{lmWinsor}
\title{
  Winsorized Regression
}
\description{
  Clip inputs and predictions to (upper, lower) or to selected quantiles
  to limit wild predictions outside the training set.
}
\usage{
  lmWinsor(formula, data, lower=NULL, upper=NULL, trim=0,
        quantileType=7, subset, weights=NULL, na.action,
        model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
        singular.ok = TRUE, contrasts = NULL, offset=NULL,
        method=c('QP', 'clip'), eps=sqrt(.Machine$double.eps),
        trace=ceiling(sqrt(nrow(data))), ...)
}
\arguments{
  \item{formula}{
    an object of class '"formula"' (or one that can be coerced to that
    class): a symbolic description of the model to be fitted.  See
    \link{lm}.  The left hand side of 'formula' must be a single vector
    in 'data', untransformed.
  }
  \item{data}{
    an optional data frame, list or environment (or object coercible by
    'as.data.frame' to a data frame) containing the variables in the
    model.  If not found in 'data', the variables are taken from
    'environment(formula)';  see \link{lm}.
  }
  \item{lower, upper}{
    Either NULL or a numeric vector or a list.

    If a numeric vector, it must have names matching columns of 'data'
    giving limits on the ranges of predictors and predictions:  If
    present, values below 'lower' will be increased to 'lower', and
    values above 'upper' will be decreased to 'upper'.  If absent, these
    limit(s) will be inferred from quantile(..., prob=c(trim, 1-trim),
    na.rm=TRUE, type=quantileType).

    If a list, its components must be numeric vectors with names, as
    just described.

    If length(trim) > 1 or either 'lower' or 'upper' is a list,
    'lmWinsor' will return a list giving alternative fits;  see 'value'
    below.
  }
  \item{trim}{
    Either a constant or a numeric vector giving the fraction (0 to 0.5)
    of observations to be considered outside the range of the data in
    determining limits not specified in 'lower' and 'upper'.  If
    length(trim) > 1 or either 'lower' or 'upper' is a list, 'lmWinsor'
    will return a list giving alternative fits;  see 'value' below.

    NOTES:

    (1) trim>0 with a singular fit may give an error.  In such cases,
    fix the singularity and retry.

    (2) trim = 0.5 should NOT be used except to check the algorithm,
    because it trims everything to the median, thereby providing zero
    leverage for estimating a regression.
  }
  \item{quantileType}{
    an integer between 1 and 9 selecting one of the nine quantile
    algorithms to be used with 'trim' to determine limits not provided
    with 'lower' and 'upper';  see   \code{\link{quantile}}.
  }
  \item{ subset }{
    an optional vector specifying a subset of observations to be used in
    the fitting process.
  }
  \item{ weights }{
    an optional vector of weights to be used in the fitting process.
    Should be 'NULL' or a numeric vector. If non-NULL, weighted least
    squares is used with weights 'weights' (that is, minimizing
    'sum(w*e*e)'); otherwise ordinary least squares is used.
  }
  \item{ na.action }{
    a function which indicates what should happen when the data contain
    'NA's.  The default is set by the 'na.action' setting of 'options',
    and is 'na.fail' if that is unset.  The factory-fresh default is
    'na.omit'.  Another possible value is 'NULL', no action.  Value
    'na.exclude' can be useful.
  }
  \item{model, x, y, qr}{
    logicals.  If 'TRUE' the corresponding components of the fit (the
    model frame, the model matrix, the response, the QR decomposition)
    are returned.
  }
  \item{ singular.ok }{
    logical. If 'FALSE' (the default in S but not in R) a singular fit
    is an error.
  }
  \item{ contrasts }{
    an optional list. See the 'contrasts.arg' of
    'model.matrix.default'.
  }
  \item{offset}{
    this can be used to specify an a priori known component to be
    included in the linear predictor during fitting.  This should be
    'NULL' or a numeric vector of length either one or equal to the
    number of cases. One or more 'offset' terms can be included in the
    formula instead or as well, and if both are specified their sum is
    used.  See 'model.offset'.
  }
  \item{method}{
    Either 'QP' or 'clip';  partial matching is allowed.  If 'clip',
    force all 'fitted.values' from an 'lm' fit to be at least
    lower[yName] and at most'upper[yName].  If 'QP', iteratively find
    the prediction farthest outside (lower, upper)[yName] and transfer
    it from the sum of squared residuals objective function to a
    constraint that keeps high 'fitted.values' from going below
    upper[yName] and low 'fitted.values' from going above lower[yName].
  }
  \item{eps}{
    small positive number used in two ways:

    \itemize{
      \item{limits}{
	'pred' is judged between 'lower' and 'upper' for 'y' as follows:
	First compute mod = mean(abs(y)).  If this is 0, let Eps = eps;
	otherwise let Eps = eps*mod.  Then pred is low if it is less
	than (lower - Eps), high if it exceeds (upper + Eps), and
	inside limits otherwise.
      }
      \item{QP}{
	To identify singularity in the quadratic program (QP) discussed
	in 'details', step 7 below, first compute the model.matrix of
	the points with interior predictions.  Then compute the QR
	decomposition of this reduced model.matix.  Then compute the
	absolute values of the diagonal elements of R.  If the smallest
	of these numbers is less than eps times the largest, terminate
	the QP with the previous parameter estimates.
      }
    }
  }
  \item{trace}{
    Print the iteration count every 'trace' iteration during 'QP'
    iterations;  see details.  'trace' = 0 means no trace.
  }
  \item{\dots}{
    additional arguments to be passed to the low level regression
    fitting functions;  see \link{lm}.
  }
}
\details{
  1.  Identify inputs and outputs via mdly <- mdlx <- formula;
  mdly[[3]] <- NULL;  mdlx[[2]] <- NULL;  xNames <- all.vars(mdlx);
  yNames <- all.vars(mdly).  Give an error if as.character(mdly[[2]]) !=
  yNames.

  2.  Do 'lower' and 'upper' contain limits for all numeric columns of
  'data?  Create limits to fill any missing.

  3.  clipData = data with all xNames clipped to (lower, upper).

  4.  fit0 <- lm(formula, clipData, subset = subset, weights = weights,
  na.action = na.action, method = method, x=x, y=y, qr=qr,
  singular.ok=singular.ok, contrasts=contrasts, offset=offset, ...)

  5.  out = a logical matrix with two columns, indicating any of
  predict(fit0) outside (lower, upper)[yName].

  6.  Add components lower and upper to fit0 and convert it to class
  c('lmWinsor', 'lm').

  7.  If((method == 'clip') || !any(out)), return(fit0).

  8.  Else, use quadratic programming (\link[quadprog]{solve.QP}) to
  minimize the 'Winsorized sum of squares of residuals' as follows:

  8.1.  First find the prediction farthest outside (lower,
  upper)[yNames].  Set temporary limits at the next closest point inside
  that point (or at the limit if that's closer).

  8.2.  Use QP to minimize the sum of squares of residuals among all
  points not outside the temporary limits while keeping the prediction
  for the exceptional point away from the interior of (lower,
  upper)[yNames].

  8.3.  Are the predictions for all points unconstrained in QP inside
  (lower, upper)[yNames]?  If yes, quit.

  8.4.  Otherwise, among the points still unconstrained, find the
  prediction farthest outside (lower, upper)[yNames].  Adjust the
  temporary limits to the next closest point inside that point (or at
  the limit if that's closer).

  8.5.  Use QP as in 8.2 but with multiple exceptional points, then
  return to step 8.3.

  9.  Modify the components of fit0 as appropriate and return the
  result.
}
\value{
  an object of class 'lmWinsor', which is either an object of class
  c('lmWinsor', 'lm') or a list of such objects.  A list is returned
  when length(trim) > 0 or is.list(lower) or is.list(upper).  The length
  of the ouput list is the max of length(trim) and the lengths of any
  lists provided for 'lower' and 'upper'.  If these lengths are not the
  same, shorter ones are replicated to match the longest one.

  An object of class c('lmWinsor', 'lm') has 'lower', 'upper', 'out',
  'message', and 'elapsed.time' components in addition to the standard
  'lm' components. The 'out' component is a logical matrix identifying
  which predictions from the initial 'lm' fit were below lower[yName]
  and above upper[yName].  If method = 'QP' and the initial fit produces
  predictions outside the limits, this object returned will also include
  a component 'coefIter' containing the model coefficients, the index
  number of the observation in 'data' transferred from the objective
  function to the constraints on that iteration, plus the sum of squared
  residuals before and after clipping the predictions and the number of
  predictions in 5 categories:  below and at the lower limit, between
  the limits, and at and above the upper limit.  The 'elapsed.time'
  component gives the run time in seconds.

  The options for 'message' are as follows:

  \item{1}{
    'Initial fit in bounds':  All predictions were between 'lower' and
    'upper' for 'y'.
  }
  \item{2}{
    'QP iterations successful':  The QP iteration described in
    'Details', step 7, terminated with all predictions either at or
    between the 'lower' and 'upper' for 'y'.
  }
  \item{3}{
    'Iteration terminated by a singular quadratic program':  The QP
    iteration described in 'Details', step 7, terminated when the
    model.matrix for the QP objective function became rank deficient.
    (Rank deficient in this case means that the smallest singular
    value is less than 'eps' times the largest.)
  }

  In addition to the coefficients, 'coefIter' also includes columns for
  'SSEraw' and 'SSEclipped', containing the residual sums of squres from
  the estimated linear model before and after clipping to the 'lower'
  and 'upper' limits for 'y', plus 'nLoOut', 'nLo.', 'nIn', 'nHi.', and
  'nHiOut', summarizing the distribtion of model predictions at each
  iteration relative to the limits.
}
\author{ Spencer Graves }
\seealso{
  \code{\link{predict.lmWinsor}}
  \code{\link{lmeWinsor}}
  \code{\link{lm}}
  \code{\link{quantile}}
  \code{\link[quadprog]{solve.QP}}
}
\examples{
# example from 'anscombe'
lm.1 <- lmWinsor(y1~x1, data=anscombe)

# no leverage to estimate the slope
lm.1.5 <- lmWinsor(y1~x1, data=anscombe, trim=0.5)

# test nonlinear optimization
lm.1.25 <- lmWinsor(y1~x1, data=anscombe, trim=0.25)

# list example
lm.1. <- lmWinsor(y1~x1, data=anscombe, trim=c(0, 0.25, .4, .5))

}
\keyword{ models }
