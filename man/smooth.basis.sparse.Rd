\name{smooth.basis.sparse}
\alias{smooth.basis.sparse}
\title{
  Construct a functional data object by smoothing data using a roughness
  penalty
}
\description{
  Makes it possible to perform smoothing with \code{smooth.basis} when the 
  data has NAs.
}
\usage{
  smooth.basis.sparse(argvals, y, fdParobj, fdnames=NULL, covariates=NULL, 
                      method="chol", dfscale=1 )
}
\arguments{
  \item{argvals}{
    a set of argument values corresponding to the observations in array
    \code{y}.  In most applications these values will be common to all curves
    and all variables, and therefore be defined as a vector or as a matrix
    with a single column.  But it is possible that these argument values
    will vary from one curve to another, and in this case \code{argvals} will
    be input as a matrix with rows corresponding to observation points and
    columns corresponding to curves.  Argument values can even vary from one
    variable to another, in which case they are input as an array with
    dimensions corresponding to observation points, curves and variables,
    respectively.  Note, however, that the number of observation points per
    curve and per variable may NOT vary.  If it does, then curves and variables
    must be smoothed individually rather than by a single call to this function.
    The default value for \code{argvals} are the integers 1 to \code{n}, where
    \code{n} is the size of the first dimension in argument \code{y}.
  }
  \item{y}{
    an set of values of curves at discrete sampling points or
    argument values. If the set is supplied as a matrix object, the rows must
    correspond to argument values and columns to replications, and it will be
    assumed that there is only one variable per observation.  If
    \code{y} is a three-dimensional array, the first dimension
    corresponds to argument values, the second to replications, and the
    third to variables within replications.  If \code{y} is a vector,
    only one replicate and variable are assumed.  If the data
    come from a single replication but multiple vectors, such as data
    on coordinates for a single space curve, then be sure to coerce
    the data into an array object by using the \code{as.array} function
    with one as the central dimension length.
  }
  \item{fdParobj}{
    a functional parameter object, a functional data object or a
    functional basis object.  In the simplest case, \code{fdParobj} may
    be a functional basis object with class "basisfd" defined
    previously by one of the "create" functions, and in this case, no
    roughness penalty is used.  No roughness penalty is also the
    consequence of supplying a functional data object with class "fd",
    in which case the basis system used for smoothing is that defining
    this object.  
    However, if the object is a functional parameter object with class
    "fdPar", then the linear differential operator object and the
    smoothing parameter in this object define the roughness penalty. For
    further details on how the roughness penalty is defined, see the help
    file for "fdPar". In general, better results can be obtained using a
    good roughness penalty than can be obtained by merely varying the
    number of basis functions in the expansion.
  }
  \item{fdnames}{
    a list of length 3 containing character vectors of names for the
    following:

    \describe{
      \item{args}{
        name for each observation or point in time at which data are
        collected for each 'rep', unit or subject.
      }
      \item{reps}{
        name for each 'rep', unit or subject.
      }
      \item{fun}{
        name for each 'fun' or (response) variable measured repeatedly
        (per 'args') for each 'rep'.
      }
    }
  }
  \item{covariates}{
    The observed values in \code{y} are assumed to be primarily
    determined by the height of the curve being estimates, but from time
    to time certain values can also be influenced by other known
    variables.  For example, multi-year sets of climate variables may be
    also determined by the presence of absence of an El Nino event, or a
    volcanic eruption.  One or more of these covariates can be supplied
    as an \code{n} by \code{p} matrix, where \code{p} is the number of
    such covariates.  When such covariates are available, the smoothing
    is called "semi-parametric."  Matrices or arrays of regression
    coefficients are then estimated that define the impacts of each of
    these covariates for each curve and each variable.
  }
  \item{method}{
    by default the function uses the usual textbook equations for computing
    the coefficients of the basis function expansions.  But, as in regression
    analysis, a price is paid in terms of rounding error for such
    computations since they involved cross-products of  basis function
    values.  Optionally, if \code{method} is set equal to the string "qr",
    the computation uses an algorithm based on the qr-decomposition which
    is more accurate, but will require substantially more computing time
    when \code{n} is large, meaning more than 500 or so.  The default
    is "chol", referring the Choleski decomposition of a symmetric positive
    definite matrix.
  }
  \item{dfscale}{
    the generalized cross-validation or "gcv" criterion that is often used
    to determine the size of the smoothing parameter involves the
    subtraction of an measue of degrees of freedom from \code{n}.  Chong
    Gu has argued that multiplying this degrees of freedom measure by
    a constant slightly greater than 1, such as 1.2, can produce better
    decisions about the level of smoothing to be used.  The default value
    is, however, 1.0.
  }
}
\value{
  an object of class \code{fd}.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{smooth.basis}}
}
