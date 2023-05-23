eval.surp <- function(evalarg, Wfdobj, nderiv=0) {
  #  Evaluates a value of a surprisal coordinate functional data object. 
  # #  Evaluates a value or a derivative of a surprisal coordinate functional  
  # #  data object. 
  #  A positive functional data object h  is <- the form
  #           h(x) <- (exp Wfdobj)(x)
  #  Note that the first two arguments may be interchanged.
  #  
  #
  #  Arguments:
  #  EVALARG ... A vector of values at which all functions are to 
  #              evaluated.
  #  WFDOBJ  ... Functional data object.  It must define a single
  #              functional data observation.
  #  NDERIV  ... The order of the derivative.  Must be 0, 1 or 2.
  #  Returns:  An array of function values corresponding to the 
  #              argument values in EVALARG
  
  #  Last modified 25 April 2025 by Jim Ramsay
  
  #  check arguments
  
  if (floor(nderiv) != nderiv) {
    stop('Third argument nderiv is not an integer.')
  }
  
  if (nderiv < 0 || nderiv > 2) {
    stop('Third argument nderiv is not 0, 1 or 2.')
  }
  
  #  Exchange the first two arguments if the first is an FD object
  #    and the second numeric
  
  if (is.numeric(Wfdobj) && fda::is.fd(evalarg)) {
    temp    <- Wfdobj
    Wfdobj  <- evalarg
    evalarg <- temp
  }
  
  #  Check the arguments
  
  if (!is.numeric(evalarg)) {
    stop('Argument EVALARG is not numeric.')
  }
  
  #  transpose EVALARG if necessary to make it a column vector
  
  evaldim <- dim(as.matrix(evalarg))
  if (evaldim[1] == 1 && evaldim[2] > 1) {  
    evalarg <- t(evalarg)  
  }
  
  #  check EVALARG
  
  if (evaldim[1] > 1 && evaldim[2] > 1) 
    stop('Argument EVALARG is not a vector.')
    
  evalarg  <- as.vector(evalarg)
  basisfd  <- Wfdobj$basis
  rangeval <- basisfd$rangeval
  evalarg[evalarg < rangeval[1]-1e-10] <- NA
  evalarg[evalarg > rangeval[2]+1e-10] <- NA
  
  #  check FDOBJ
    
  if (!fda::is.fd(Wfdobj)) 
        stop('Argument FD is not a functional data object.')
  
  #  compute Zmat, a M by M-1 orthonormal matrix
  
  Bmat <- Wfdobj$coefs
  M    <- dim(Bmat)[2] + 1
  if (M == 2) {
    root2 <- sqrt(2)
    Zmat <- matrix(1/c(root2,-root2),2,1)
  } else {
    Zmat <- zerobasis(M)
  }
  
  #  Set up coefficient array for FD
  
  BZmat    <- Bmat %*% t(Zmat)
  Xmat     <- fda::eval.basis(evalarg, basisfd) %*% BZmat
  MXmat    <- M^Xmat
  sum0     <- rowSums(MXmat)
  SumMXmat <- matrix(rep(sum0,each=M), ncol=M, byrow=TRUE)
  
  #  Case where EVALARG is a vector of values to be used for all curves
  #  NB:  Resist the temptation of use /log(M) anywhere else.  
  #  This code sets up Smat as defined with log basis M, and no further
  #  definition of log basis is required.
  
  if (nderiv == 0) {
    Smat   <- -Xmat + log(SumMXmat)/log(M)
    return(Smat)
  }
  
  # Note:  derivative values computed using Pmat rather than Bmat
  # This needs correcting
  
  #  First derivative:
  
  if (nderiv == 1) {
    Pmat    <- MXmat/SumMXmat
    DXmat   <- fda::eval.basis(evalarg, basisfd, 1) %*% BZmat
    matSum  <- rowSums(Pmat*DXmat)
    Rmat    <- matrix(rep(matSum,each=M), ncol=M, byrow=TRUE)
    DSmat   <- -DXmat + Rmat
    return(DSmat)
  }
  
  #  Second derivative
  
  if (nderiv == 2) {
    Pmat    <- MXmat/SumMXmat
    DXmat   <- fda::eval.basis(evalarg, basisfd, 1) %*% BZmat
    D2Xmat  <- fda::eval.basis(evalarg, basisfd, 2) %*% BZmat
    matSum  <- rowSums(Pmat*DXmat)
    Rmat    <- matrix(rep(matSum,each=M), ncol=M, byrow=TRUE)
    matSum2 <- rowSums((Pmat*(D2Xmat + DXmat^2) - Rmat*DXmat))
    D2Smat  <- -D2Xmat + 
      matrix(rep(matSum2,each=M), ncol=M, byrow=TRUE)
    return(D2Smat)
  }
}

