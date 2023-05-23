smooth.basis.sparse <- function(argvals, y, fdParobj, fdnames=NULL, covariates=NULL, 
                                method="chol", dfscale=1 ){
  
  #  Arguments:
  # ARGVALS  A set of N argument values, set by default to equally spaced
  #             on the unit interval (0,1).
  # Y        an array containing values of curves
  #             If the array is a matrix, rows must correspond to argument
  #             values and columns to replications, and it will be assumed
  #             that there is only one variable per observation.
  #             If Y is a three-dimensional array, the first dimension
  #             corresponds to argument values, the second to replications,
  #             and the third to variables within replications.
  #             If Y is a vector, only one replicate and variable are assumed.
  # FDPAROBJ A functional parameter or fdPar object.  This object
  #             contains the specifications for the functional data
  #             object to be estimated by smoothing the data.  See
  #             comment lines in function fdPar for details.
  #             This argument may also be either a FD object, or a
  #             BASIS object.  In this case, the smoothing parameter
  #             LAMBDA is set to 0.
  # FDNAMES  A cell of length 3 with names for
  #             1. argument domain, such as "Time"
  #             2. replications or cases
  #             3. the function.
  # COVARIATES  A N by Q matrix Z of covariate values used to augment
  #             the smoothing function, where N is the number of
  #             data values to be smoothed and Q is the number of
  #             covariates.  The process of augmenting a smoothing
  #             function in this way is often called "semi-parametric
  #             regression".  The default is the null object NULL.
  # METHOD      The method for computing coefficients.  The usual method
  #             computes cross-product matrices of the basis value matrix,
  #             adds the roughness penalty, and uses the Choleski
  #             decomposition of this to compute coefficients, analogous
  #             to using the normal equations in least squares fitting.
  #             But this approach, while fast, contributes unnecessary
  #             rounding error, and the qr decomposition of the augmented
  #             basis matrix is prefererable.  But nothing comes for free,
  #             and the computational overhead of the qr approach can be a
  #             serious problem for large problems (n of 1000 or more).
  #             For this reason, the default is "method" = "chol", but if
  #             'method' == 'qr', the qr decomposition is used.
  # DFFACTOR A multiplier of df in GCV, set to one by default
  #
  # Returns a list containing:
  #   FDOBJ   an object of class fd containing coefficients.
  #   DF      a degrees of freedom measure.
  #   GCV     a measure of lack of fit discounted for df.
  #              If the function is univariate, GCV is a vector
  #              containing the error  sum of squares for each
  #              function, and if the function is multivariate,
  #              GCV is a NVAR by NCURVES matrix.
  #   COEF    the coefficient matrix for the basis function
  #                expansion of the smoothing function
  #   SSE     the error sums of squares.
  #              SSE is a vector or matrix of the same size as
  #              GCV.
  #   PENMAT  the penalty matrix.
  #   Y2CMAP  the matrix mapping the data to the coefficients.
  
  if (is.fdPar(fdParobj)) {
    basisobj = fdParobj$fd$basis
  } else {
    if (is.fd(fdParobj)) {
      basisobj = fdParobj$basis
    } else {
      if (is.basis(fdParobj)) {
        basisobj = fdParobj
      } else {
        stop("fdParobj is not a fdPar, fd, or a basis object.")
      }
    }
  }
  if(length(dim(y)) == 2){
	coefs = matrix(0, nrow = basisobj$nbasis, ncol = dim(y)[2])
	for(i in 1:dim(y)[2]){
		curve = y[,i]
		curve.smooth = smooth.basis(argvals[!is.na(curve)],curve[!is.na(curve)],
                                basisobj, covariates, method)
		coefs[,i] = curve.smooth$fd$coefs
	}
  } else if(length(dim(y)) == 3){
 	coefs = array(0, c(basisobj$nbasis,dim(y)[2:3]))
	for(i in 1:dim(y)[2]){
		for(j in 1:dim(y)[3]){
			curve = y[,i,j]
			curve.smooth = smooth.basis(argvals[!is.na(curve)],curve[!is.na(curve)],
                                basisobj, covariates, method)
			coefs[,i,j] = curve.smooth$fd$coefs
		}
	} 
  }
  datafd = fd(coefs,basisobj, fdnames)
  return(datafd)
}