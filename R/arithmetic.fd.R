'+.fd' <- function(fd1, fd2)
{
	#  Arguments:
	#  FDA1  ...  A functional data object
	#  FDA2  ...  A functional data object
	#  Returns:
	#  FDASUM  ...  A functional data object that is FDA1 plus FDA2
	#  Last modified:  4 July 2001
	
    if (!(inherits(fd1, "fd") & inherits(fd2, "fd"))) stop(
      'Both arguments to + are not functional data objects.')
    coef1 <- getcoef(fd1)
    coef2 <- getcoef(fd2)
    if (any(dim(coef1) != dim(coef2))) stop(
      'Coefficient arrays are not of same dimensions for +.')
    fdsum <- create.fd(coef1+coef2, getbasis(fd1), fd1$fdnames)
    return(fdsum)
}

#  ----------------------------------------------------------------

'-.fd' <- function(fd1, fd2)
{
	#  Arguments:
	#  FDA1  ...  A functional data object
	#  FDA2  ...  A functional data object
	#  Returns:
	#  FDADIF  ...  A functional data object that is FDA1 minus FDA2
	#  Last modified:  4 July 2001
	
	if (!(inherits(fd1, "fd") & inherits(fd2, "fd"))) stop(
      'Both arguments to + are not functional data objects.')
    coef1 <- getcoef(fd1)
    coef2 <- getcoef(fd2)
    if (any(dim(coef1) != dim(coef2))) stop(
      'Coefficient arrays are not of same dimensions for -.')
    fddif <- create.fd(coef1-coef2, getbasis(fd1), fd1$fdnames)
    return(fddif)
}

#  ----------------------------------------------------------------

'*.fd' <- function(fd1, fd2)
{
 	#  Arguments:
	#  FDA1  ...  Either a functional data object or a number
	#  FDA2  ...  Either a functional data object or a number
	#  At least oneof FDA1 and FDA2 must be a functional data object.
	#  Returns:
	#  FDAPROD  ...  A functional data object that is FDA1 times FDA2
	#  Last modified:  4 July 2001
	
   if ((!(inherits(fd1, "fd") | inherits(fd2, "fd")))) stop(
      'Neither argument for * is a functional data object.')
    if (inherits(fd1, "fd") & inherits(fd2, "fd") ) {
    coef1     <- getcoef(fd1)
    coefd1    <- dim(coef1)
    coef2     <- getcoef(fd2)
    coefd2    <- dim(coef2)
    if (length(coefd1) != length(coefd2)) stop(
      'Number of dimensions of coefficient arrays do not match.')
    if (any(coefd1 != coefd2)) stop(
      'Dimensions of coefficient arrays do not match.')
    basisfd1  <- getbasis(fd1)
    basisfd2  <- getbasis(fd2)
    nbasis1   <- basisfd1$nbasis
    nbasis2   <- basisfd2$nbasis
    rangeval1 <- basisfd1$rangeval
    rangeval2 <- basisfd2$rangeval
    if (any(rangeval1 != rangeval2)) stop(
      'The ranges of the arguments are not equal.')
    neval     <- max(10*max(nbasis1,nbasis2) + 1, 101)
    evalarg   <- seq(rangeval1[1],rangeval2[2], length=neval)
    fdarray1  <- eval.fd(evalarg, fd1)
    fdarray2  <- eval.fd(evalarg, fd2)
    fdarray   <- fdarray1*fdarray2
    if (nbasis1 > nbasis2) {
      basisfd  <- basisfd1
      coefprod <- project.basis(fdarray, evalarg, basisfd)
    } else {
      basisfd  <- basisfd2
      coefprod <- project.basis(fdarray, evalarg, basisfd)
    }
    fdnames1 <- fd1$fdnames
    fdnames2 <- fd2$fdnames
    fdnames  <- fdnames1
    fdnames[[3]] <- paste(fdnames1[[3]],'*',fdnames2[[3]])
    fdprod   <- create.fd(coefprod, basisfd, fdnames)
    return(fdprod)
  } else {
    if ((!(is.numeric(fd1) | is.numeric(fd2)))) stop(
      'Neither argument for * is numeric.')
    if (is.numeric(fd1)) {
      fac <- fd1
      fd  <- fd2
    } else {
      fac <- fd2
      fd  <- fd1
    }
    coef     <- getcoef(fd)
    coefd    <- dim(coef)
    basisfd  <- getbasis(fd)
    nbasis   <- basisfd$nbasis
    rangeval <- basisfd$rangeval
    neval    <- max(10*nbasis + 1,101)
    evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
    fdarray  <- fac*eval.fd(evalarg, fd)
    coefprod <- project.basis(fdarray, evalarg, basisfd)
    fdnames  <- fd$fdnames
    fdnames[[3]] <- paste(as.character(fac),'*',fdnames[[3]])
    fdprod   <- create.fd(coefprod, basisfd, fdnames)
    return(fdprod)
 }
}

#  ----------------------------------------------------------------

'/.fd' <- function(fd1, fd2)
{
  	#  Arguments:
	#  FDA1  ...  Either a functional data object or a number
	#  FDA2  ...  Either a functional data object or a number
	#  At least oneof FDA1 and FDA2 must be a functional data object.
	#  Returns:
	#  FDAPROD  ...  A functional data object that is FDA1 divided by FDA2
	#  Last modified:  4 July 2001
	
   if ((!(inherits(fd1, "fd") | inherits(fd2, "fd")))) stop(
      'Neither argument for * is a functional data object.')
    if (inherits(fd1, "fd") & inherits(fd2, "fd")) {
    coef1     <- getcoef(fd1)
    coefd1    <- dim(coef1)
    coef2     <- getcoef(fd2)
    coefd2    <- dim(coef2)
    if (length(coefd1) != length(coefd2)) stop(
      'Number of dimensions of coefficient arrays do not match.')
    if (any(coefd1 != coefd2)) stop(
      'Dimensions of coefficient arrays do not match.')
    basisfd1  <- getbasis(fd1)
    basisfd2  <- getbasis(fd2)
    nbasis1   <- basisfd1$nbasis
    nbasis2   <- basisfd2$nbasis
    rangeval1 <- basisfd1$rangeval
    rangeval2 <- basisfd2$rangeval
    if (any(rangeval1 != rangeval2)) stop(
      'The ranges of the arguments are not equal.')
    neval     <- max(10*max(nbasis1,nbasis2) + 1, 101)
    evalarg   <- seq(rangeval1[1],rangeval2[2], length=neval)
    fdarray1  <- eval.fd(evalarg, fd1)
    fdarray2  <- eval.fd(evalarg, fd2)
    fdarray   <- fdarray1/fdarray2
    if (nbasis1 > nbasis2) {
      basisfd  <- basisfd1
      coefquot <- project.basis(fdarray, evalarg, basisfd)
    } else {
      basisfd  <- basisfd2
      coefquot <- project.basis(fdarray, evalarg, basisfd)
    }
    fdnames1 <- fd1$fdnames
    fdnames2 <- fd2$fdnames
    fdnames  <- fdnames1
    fdnames[[3]] <- paste(fdnames1[[3]],'*',fdnames2[[3]])
    fdquot   <- create.fd(coefquot, basisfd, fdnames)
    return(fdquot)
  } else {
    if ((!(is.numeric(fd1) | is.numeric(fd2)))) stop(
      'Neither argument for * is numeric.')
    if (is.numeric(fd1)) {
      fac <- fd1
      fd  <- fd2
    } else {
      fac <- fd2
      fd  <- fd1
    }
    coef     <- getcoef(fd)
    coefd    <- dim(coef)
    basisfd  <- getbasis(fd)
    nbasis   <- basisfd$nbasis
    rangeval <- basisfd$rangeval
    neval    <- max(10*nbasis + 1, 101)
    evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
    fdnames  <- fd$fdnames
    if (is.numeric(fd1)) {
      fdarray  <- fac/eval.fd(evalarg, fd)
      fdnames[[3]] <- paste(as.character(fac),'/',fdnames[[3]])
    } else {
      fdarray  <- eval.fd(evalarg, fd)/fac
      fdnames[[3]] <- paste(fdnames[[3]],'/',as.character(fac))
    }
    coefquot <- project.basis(fdarray, evalarg, basisfd)
    fdquot   <- create.fd(coefquot, basisfd, fdnames)
    return(fdquot)
 }
}

#  ----------------------------------------------------------------

'^.fd' <- function(fd, power)
{
   	#  Arguments:
	#  FD     ...  A functional data object
	#  POWER  ...  An exponent
	#  Returns:
	#  FDAPOWR  ...  A functional data object that is FD to the power POWER
	#  Last modified:  4 July 2001
	
  if ((!(inherits(fd, "fd")))) stop(
      'First argument for ^ is not a functional data object.')
  if ((!(is.numeric(power)))) stop(
      'Second argument for ^ is not numeric.')
  coef     <- getcoef(fd)
  coefd    <- dim(coef)
  basisfd  <- getbasis(fd)
  nbasis   <- basisfd$nbasis
  rangeval <- basisfd$rangeval
  neval    <- max(10*nbasis + 1,101)
  evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
  fdnames  <- fd$fdnames
  fdarray  <- eval.fd(evalarg, fd)^power
  coefpowr <- project.basis(fdarray, evalarg, basisfd)
  fdpowr   <- create.fd(coefpowr, basisfd, fdnames)
  return(fdpowr)
}


#  ----------------------------------------------------------------

'sqrt.fd' <- function(fd)
{
   	#  Arguments:
	#  FD     ...  A functional data object
 	#  Returns:
	#  FDASQRT  ...  A functional data object that is the square root of FD
	#  Last modified:  4 July 2001
	
  if ((!(inherits(fd, "fd")))) stop(
      'First argument is not a functional data object.')
  coef     <- getcoef(fd)
  coefd    <- dim(coef)
  basisfd  <- getbasis(fd)
  nbasis   <- basisfd$nbasis
  rangeval <- basisfd$rangeval
  neval    <- max(10*nbasis + 1,101)
  evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
  fdnames  <- fd$fdnames
  fdarray  <- sqrt(eval.fd(evalarg, fd))
  coefsqrt <- project.basis(fdarray, evalarg, basisfd)
  fdpowr   <- create.fd(coefsqrt, basisfd, fdnames)
  return(fdsqrt)
}
