data2fd <- function(y, argvals = seq(0, 1, len = n), basisobj,
                    fdnames = defaultnames, 
                    argnames = c("time", "reps", "values") )
{
#  DATA2FD Converts an array Y of function values plus an array
#    ARGVALS of argument values into a functional data object.
#
#  A functional data object is a sample of one or more functions, called
#    functional data observations.
#  A functional data observation consists of one or more
#    curves, each curve corresponding to a variable.
#  For example, a functional data object can be a sample of size 35 of
#    temperature functions, one for each of 35 Canadian weather stations.
#    In this case, each observations consists of a single temperature
#    function.
#  Or, for example, a functional data object can be a sample of size 35
#    of temperature and precipitation functional observations.  In this case
#    each observations consists of two curves, one for the temperature
#    and one for the precipitation variable.
#  All functional objects have a one-dimensional argument.  In the above
#    examples, this argument is time measured in months or days.
#
#  The data values used by DATA2FD to make the functional data object FDOBJ
#    are in array Y.  If each observation consists of a single curve,
#    Y is two-dimensional If each observations has multiple variables,
#    Y will be three-dimensional.  See below for further details.
#
#  DATA2FD assumes that each observation is evaluated at a set of argument
#    specified in ARGVALS. These argument values may be common to all curves,
#    in which case ARGVALS is a one-dimensional vector.  If argument values
#    vary from observation to observation, ARGVALS will be two dimensional.
#    See below for further details.
#
#  Arguments for this function are as follows.  The first three are necessary
#    and the fourth is optional.
#
#  Y ... (necessary)
#  The array Y stores curve values used to create functional data object FDOBJ.
#  Y can have one, two, or three dimensions according to whether whether
#    the sample size, the number of variables in each observation.  Its
#    dimensions are:
#     1.  argument values  ...  size = no. argument values in ARGVAL
#     2.  replications     ...  size = sample size
#     3.  variables        ...  size = no. variables per observation
#  If Y is a one-way array, either as a vector or a matrix with one column,
#     it's single non-trivial dimension = no. argument values.  If Y
#     is two-dimensional, each observation is assumed to have one variable.
#     If Y is three-dimensional, each observation is assumed to have
#     multiple variables.  Note:  a single multivariate observation must
#     be an array Y with three dimensions, the middle of which is of length 1.
#  The values in Y may be missing, indicated by NaN.  The presence of missing
#     values will slow down the computation somewhat since each observation
#     must then be processed separately.
#  Example:  For monthly temperature data for 35 weather stations,
#     Y will be 12 by 35.  For both temperature and precipitation observations,
#     Y will be 12 by 35 by 2.  For temperature/precipitation data at Montreal
#     only, Y will be 12 by 1 by 2.
#  This argument is necessary, and there is no default value.
#
#  ARGVALS  ... (necessary)
#  A set of argument values.  In most situations, these will be common to all
#    observations, and ARGVALS will be a one-dimensional vector, with one
#    element per observation.  These values need not be increasing.
#    In the weather station example for monthly data, ARGVALS is
#    a vector of length 12, with values 0.5, 1.5,..., 11.5.
#    However, ARGVALS may also be a matrix if the argument values vary from
#    observation to observation.  In this case, the second dimension is
#    the same as that of Y.  If the number of arguments also varies from
#    from observation to observation, the second dimension of ARGVALS should
#    equal the the largest number of argument values, and the elements in
#    each show should be padded out with NaN.
#    Argument values falling outside of the range specified in the
#    BASISOBJ and their corresponding values in Y will not be used,
#    but if this happens, a warning message is displayed.
#  This argument is necessary, and there is no default value.
#
#  BASISOBJ  ...  (necessary)
#    A functional data basis object created by function CREATE.BASIS.FD
#    or one of its specialized version, such as CREATE.BSPLINE.BASIS or
#    CREATE.FOURIER.BASIS.  The functional data basis object specifies
#    a basis type (eg. 'fourier' or 'bspline'), a range or argument values,
#    the number of basis functions, and fixed parameters determining these
#    basis functions (eg. period for 'fourier' bases or knots for 'bspline'
#    bases.
#    In most applications, BASISOBJ will be supplied.  If BASISOBJ is supplied,
#    the next three arguments are ignored.
#    If BASISOBJ is an essential argument, and there no default value.  But
#    see function MAKE.BASIS for a simplified technique for defining this
#    basis.  For example, function call
#         MAKE.BASIS([0,12], 7, 1)
#    could be used for the monthly temperature/precipitation data to define
#    a 'fourier' basis over an interval of 12 months containing 7 basis
#    functions (the 3rd argument species the basis to be periodic.)
#    This argument is necessary, and there is no default value.
#    Earlier releases of DATA2FD supplied additional arguments for constructing
#    a default basis, but these have been eliminated in favor of using new
#    function MAKE.BASIS.
#
#  FDNAMES  ... (optional)
#    A list of length 3 with members containing
#               1. a vector of names for the argument values
#               2. a vector of names for the replications or cases
#               3. a name for the function, or a vector of names if there
#                  are multiple functions.
#    For example, for the monthly temperature/precipitation data,
#    fdnames{1} = 'Month'
#    fdnames{2} = 'Station'
#    fdnames{3} = 'Variable'
#    By default, the string 'time', 'reps' and 'values' are used.
#
#  ARGNAMES ...
#    A vector  of type "character" of length 3 containing
#               1. the name of the argument, e.g. "time" or "age"
#               2. a description of the cases, e.g. "weather stations"
#               3. the name for the function(s), e.g."temperature" or "position"
#               These names are used as names for the members of list FDNAMES.
#
#  DATA2FD Returns the object FDOBJ of functional data class containing
#    coefficients for the expansion and the functional data basis object
#    basisobj.
#
#  DATA2FD is intended for more casual use not requiring a great deal of
#    control over the smoothness of the functional data object.  It uses
#    function PROJECT.BASIS to compute the functional data object. Indeed,
#    in the simplest and most common situation, DATA2FD consists of
#             coef  = project.basis(y, argvals, basisobj)
#             fdobj = fd(coef,basisobj,fdnames)
#    However, for more advanced applications requiring more smoothing
#    control than is possible by setting the number of basis functions in
#    basisobj, function SMOOTH.BASIS should be used.  Or, alternatively,
#    DATA2FD may first be used with a generous number of basis functions,
#    followed by smoothing using function SMOOTH.

#  Last modified:  26 October 2005

#
#  set up default fdnames, using dimnames of Y if there are any.
#

  defaultnames      <- vector("list",3)

#
#  check basisobj argument.
#
  if(!(inherits(basisobj, "basisfd"))) stop(
      "BASISOBJ is not a functional data basis.")
  nbasis <- basisobj$nbasis

#
#  Make Y an array, and determine its dimensions
#
  if (is.array(y) == FALSE) y <- as.array(y)
  yd   <- dim(y)
  ndim <- length(yd)
  if (ndim == 1) {
	 y <- as.matrix(y)
	 yd <- dim(y)
	 ndim <- length(yd)
  }
  if (ndim > 3) stop(
      "Too many dimensions for argument Y.")
  #  Determine the maximum number of argument values, number of replicates, and
  #    number of variables
  n  <- yd[1]
  if (n == 1) stop(
      "Only one argument value not allowed.")
  if (ndim > 1) nrep <- yd[2] else nrep <- 1
  if (ndim > 2) nvar <- yd[3] else nvar <- 1

#  set up defaultnames to be used for fdnames slot

  if (is.null(defaultnames[[3]])) defaultnames[[3]] <- as.character(1:nvar)
  names(defaultnames)[1:length(argnames)] <- argnames

#  Make ARGVALS an array and check for compatibility with Y

  if(is.array(argvals) == FALSE) argvals <- as.array(argvals)
  argd  <- dim(argvals)
  nargd <- length(argd)
  if (nargd > 2) stop(
     "ARGVALS has too many dimensions.")
  if (argd[1] != n) stop(
    "Number of argument values incompatible with number of data.")
  if (nargd==2 && argd[2] != nrep) stop(
    paste("Matrix argvals must have same number of columns\n",
          "as the number of replicates."))
#  Issue a warning of arguments are outside of the in the basisobj.
  rangeval <- basisobj$rangeval
  temp <- c(argvals)
  temp <- temp[!is.na(temp)]
  if (min(temp) < rangeval[1] | max(temp) > rangeval[2]) {
    warning(c("Some arguments values are outside of the range in basisobj,\n",
          " and some data values will not be used."))
    if (nargd == 1) {
       index <- argvals < rangeval[1] | argvals > rangeval[2]
       argvals[index] <- NA
    } else {
       for (irep in 1:nrep) {
          index <- argvals[,irep] < rangeval[1] | argvals[,irep] > rangeval[2]
          argvals[index,irep] <- NaN
       }
    }
  }
#
# Process data differently according to whether:
#
#  First case:  ARGVALS is a vector, and no missing values in y.
#  Second case: ARGVALS is a vector, but missing values in y; in this case
#      ARGVALS is turned into a matrix, but a lot of the work
#      in PROJECT.BASIS is not repeated.
#  Third case:  ARGVALS is a matrix.
#
  if(nargd == 1) {
    if (sum(is.na(argvals)) == 0 && sum(is.na(y)) == 0) {
#  First case:  no missing values, ARGVALS a vector
#    In this case all of the work is done by function PROJECT.BASIS
      if (nbasis <= n) {
        coef <- project.basis(y, argvals, basisobj)
      } else {
        coef <- project.basis(y, argvals, basisobj, TRUE)
      }
    } else {
#  Second case: ARGVALS a vector, but missing data present
      coefd    <- yd
      coefd[1] <- nbasis
      coef     <- array(dim = coefd)
#
# set up penalty and basis matrices
#
      index    <- !is.na(argvals)
      basismat <- eval.basis(argvals, basisobj)
      penmat   <- getbasispenalty(basisobj)
      # add a small amount to diagonal of penalty to ensure conditioning
      penmat   <- penmat + 1e-10 * max(penmat) * diag(dim(penmat)[1])
      # add a small penalty to deal with underdetermination by data
      lambda   <- 0.0001 * sum(basismat[index,]^2)/sum(diag(penmat))
      penmat <- lambda * penmat
      if(length(coefd) == 2) {
        #  Univariate functions
        for(j in (1:nrep)) {
          yy    <- y[,j]
          index <- !is.na(yy) & !is.na(argvals)
          if (length(yy[index]) < 2) stop(
              paste("Less than 2 data values available for curve",
                    j,"."))
          Cmat  <- crossprod( basismat[index,] ) + penmat
          Dmat  <- crossprod( basismat[index,],yy[index])
          coef[, j] <- symsolve( Cmat, Dmat )
        }
      } else {
        #  Multivariate functions
        for(j in (1:nrep))  for(k in (1:nvar)) {
          yy <- y[, j, k]
          index <- !is.na(yy) & !is.na(argvals)
          if (length(yy[index]) < 2) stop(
              paste("Less than 2 data values available for curve",
                    j," and variable",k,"."))
          Cmat <- crossprod( basismat[index,] ) + penmat
          Dmat <- crossprod( basismat[index,],yy[index] )
          coef[, j, k] <- symsolve( Cmat, Dmat )
        }
      }
    }
  } else {
#
#  now deal with the case where argvals is a matrix (no missing data)
#
    coefd    <- yd
    coefd[1] <- nbasis
    coef     <- array(dim = coefd)
    argv     <- c(argvals)
    index    <- !is.na(argv)
    argv     <- unique(argv[index])
    basismat <- eval.basis(argv, basisobj)
    penmat   <- getbasispenalty(basisobj)
    penmat   <- penmat + 1e-10 * max(penmat) * diag(dim(penmat)[1])
    # add a small penalty to deal with underdetermination by data
    lambda1  <- 0.0001 * sum(basismat^2)/sum(diag(penmat))
    penmat1  <- lambda1 * penmat
#
    if(length(coefd) == 2) {
      #  Univariate functions
      for(j in (1:nrep)) {
        yy    <- y[, j]
        argv  <- argvals[, j]
        index <- (!is.na(yy) & !is.na(argv))
        if (length(yy[index]) < 2) stop(
              paste("Less than 2 data values available for curve",
                    j,"."))
        coef[, j] <-
            project.basis(yy[index], argv[index], basisobj, TRUE)
      }
    } else {
      #  Multivariate functions
      for(j in (1:nrep))  for(k in (1:nvar)) {
        yy <- y[, j, k]
        argv <- argvals[, j]
        index <- (!is.na(yy) & !is.na(argv))
        if (length(yy[index]) < 2) stop(
              paste("Less than 2 data values available for curve",
                    j," and variable",k,"."))
        coef[, j, k] <-
            project.basis(yy[index], argv[index], basisobj, TRUE)
      }
    }
  }
  #
  #  Now that coefficient array has been computed, create functional data object
  #
  fd <- fd(coef, basisobj, fdnames = fdnames)

  fd
}

