#  setClass for "fd"

# setClass("fd",    representation(coefs    = "array",
#                                  basis    = "basisfd",
#                                  fdnames  = "list"))

#  Generator function of class fd

fd <- function(coef=matrix(0,2,1), basisobj=basisfd(), fdnames=defaultnames)
{
  #  This function creates a functional data object.
  #    A functional data object consists of a basis for expanding a functional
  #    observation and a set of coefficients defining this expansion.
  #    The basis is contained in a "basisfd" object; that is, a realization
  #    of the "basisfd" class.

  #  Arguments
  #  COEF ... An array containing coefficient values for the expansion of each
  #             set of function values in terms of a set of basis functions.
  #           If COEF is a three-way array, then the first dimension
  #             corresponds to basis functions, the second to replications,
  #             and the third to variables.
  #           If COEF is a matrix, it is assumed that there is only
  #             one variable per replication, and then
  #                 rows    correspond to basis functions
  #                 columns correspond to replications
  #           If COEF is a vector, it is assumed that there is only one
  #             replication and one variable.
  #  BASISOBJ ... a functional data basis object
  #  FDNAMES  ... The analogue of the dimnames attribute of an array, this is
  #               a list of length 3 with members containing:
  #               1. a character vector of names for the argument values
  #               2. a character vector of names for the replications or cases
  #               3. a character vector of names for the functions
  #               Each of these vectors can have a name referring to the modality
  #                 of the data.  An example would be "time", "reps", "values"

  #  Returns:
  #  FD ... a functional data object

  #  last modified 30 January 2007

    #  check basisobj

    if (!(inherits(basisobj, "basisfd"))) stop(
        "Argument basis must be of basis class")

    type <- basisobj$type

    #  check COEF and get its dimensions

    if (!is.numeric(coef)) stop("coef must be numerical vector or matrix")
    else if (is.vector(coef)) {
            coef  <- as.matrix(coef)
            if (type %in% c("const", "constant")) coef <- t(coef)
# changed from =="constant" to %in% c(...) 
# 2007.01.30 Spencer Graves            
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
    else if (is.matrix(coef)) {
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
    else if (is.array(coef)) {
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
    else stop("argument coef is not correct")

    if (ndim > 3) stop(
        "First argument not of dimension 1, 2 or 3")

    if (coefd[1] != basisobj$nbasis) stop(
        "Number of coefficients does not match number of basis functions.")

    #  setup number of replicates and number of variables

    if (ndim > 1) nrep <- coefd[2] else nrep <- 1
    if (ndim > 2) nvar <- coefd[3] else nvar <- 1

    #  set up default fdnames

    if (ndim == 1) defaultnames <- list("time", "reps", "values")
    if (ndim == 2) defaultnames <- list("time",
                                        paste("reps",as.character(1:nrep)),
                                        "values")
    if (ndim == 3) defaultnames <- list("time",
                                        paste("reps",as.character(1:nrep)),
                                        paste("values",as.character(1:nvar)) )

    names(defaultnames) <- c("args", "reps", "funs")

#  S4 definition
#	fdobj <- new("fd", coefs=coef, basis=basisobj, fdnames=fdnames)

#  S3 definition

	fdobj <- list(coefs=coef, basis=basisobj, fdnames=fdnames)
	oldClass(fdobj) <- "fd"
	
	fdobj
}

#  "print" method for "fd"

print.fd <- function(x, ...)
{
	
	cat("Functional data x:\n")
	cat(" Dimensions of the data:\n")
	cat("", paste(names(x$fdnames), collapse=", "), "\n")
	
	print.basisfd(x$basis)
	
}


#  plus method for "fd"

"+.fd" <- function(fd1, fd2)
{
	#  Arguments:
	#  FD1  ...  A functional data object
	#  FD2  ...  A functional data object
	#  Returns:
	#  FDDIF  ...  A functional data object that is FD1 plus FD2
	#  Last modified:  25 October 2005

  coef1 <- fd1$coefs
  coef2 <- fd2$coefs
  if (any(dim(coef1) != dim(coef2))) stop("Coefficient arrays are not of same dimensions for +.fd")
  fdsum <- fd(coef1+coef2, fd1$basis, fd1$fdnames)
  fdsum
}

#  minus method for "fd"

"-.fd" <- function(fd1, fd2)
{
	#  Arguments:
	#  FD1  ...  A functional data object
	#  FD2  ...  A functional data object
	#  Returns:
	#  FDDIF  ...  A functional data object that is FD1 plus FD2
	#  Last modified:  17 September 2005

  coef1 <- fd1$coefs
  coef2 <- fd2$coefs
  if (any(dim(coef1) != dim(coef2))) stop("Coefficient arrays are not of same dimensions for -.fd")
  fddif <- fd(coef1-coef2, fd1$basis, fd1$fdnames)
  fddif
}

#  point-wise product method for "fd"
"*.fd" <- function(fd1, fd2)
{
	#  Arguments:
	#  FD1  ...  Either a functional data object or a number
	#  FD2  ...  Either a functional data object or a number
	#  At least one of FD1 and FD2 must be a functional data object.
	#  Returns:
	#  FDAPROD  ...  A functional data object that is FD1 times FD2
	#  Last modified:  17 September 2005
	
   if ((!(inherits(fd1, "fd") | inherits(fd2, "fd")))) stop(
      "Neither argument for * is a functional data object.")
    if (inherits(fd1, "fd") & inherits(fd2, "fd") ) {
    	coef1     <- fd1$coefs
    	coefd1    <- dim(coef1)
    	coef2     <- fd2$coefs
    	coefd2    <- dim(coef2)
    	if (length(coefd1) != length(coefd2)) stop(
      		"Number of dimensions of coefficient arrays do not match.")
    	if (any(coefd1 != coefd2)) stop(
      		"Dimensions of coefficient arrays do not match.")
    	basisfd1  <- fd1$basis
    	basisfd2  <- fd2$basis
   		nbasis1   <- basisfd1$nbasis
    	nbasis2   <- basisfd2$nbasis
    	rangeval1 <- basisfd1$rangeval
    	rangeval2 <- basisfd2$rangeval
    	if (any(rangeval1 != rangeval2)) stop(
      		"The ranges of the arguments are not equal.")
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
    	fdnames[[3]] <- paste(fdnames1[[3]],"*",fdnames2[[3]])
    	fdprod   <- fd(coefprod, basisfd, fdnames)
    	return(fdprod)
  	} else {
    	if ((!(is.numeric(fd1) | is.numeric(fd2)))) stop(
      		"Neither argument for * is numeric.")
    	if (is.numeric(fd1)) {
      		fac <- fd1
      		fd  <- fd2
    	} else {
      		fac <- fd2
      		fd  <- fd1
    	}
    	coef     <- fd$coefs
    	coefd    <- dim(coef)
    	basisfd  <- fd$basis
    	nbasis   <- basisfd$nbasis
    	rangeval <- basisfd$rangeval
    	neval    <- max(10*nbasis + 1,101)
    	evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
    	fdarray  <- fac*eval.fd(evalarg, fd)
    	coefprod <- project.basis(fdarray, evalarg, basisfd)
    	fdnames  <- fd$fdnames
    	fdnames[[3]] <- paste(as.character(fac),"*",fdnames[[3]])
    	fdprod   <- fd(coefprod, basisfd, fdnames)
    	return(fdprod)
 	}
}

#  point-wise quotient method for "fd"

"/.fd" <- function(fd1, fd2)
{
	#  Arguments:
	#  FD1  ...  Either a functional data object or a number
	#  FD2  ...  Either a functional data object or a number
	#  At least one of FD1 and FD2 must be a functional data object.
	#  Returns:
	#  FDAQUOT  ...  A functional data object that is FD1 times FD2
	#  Last modified:  17 September 2005
	
   if ((!(inherits(fd1, "fd") | inherits(fd2, "fd")))) stop(
      "Neither argument for * is a functional data object.")
    if (inherits(fd1, "fd") & inherits(fd2, "fd") ) {
	    #  both arguments are functional data objects
    	coef1     <- fd1$coefs
    	coefd1    <- dim(coef1)
    	coef2     <- fd2$coefs
    	coefd2    <- dim(coef2)
    	if (length(coefd1) != length(coefd2)) stop(
      		"Number of dimensions of coefficient arrays do not match.")
    	if (any(coefd1 != coefd2)) stop(
      		"Dimensions of coefficient arrays do not match.")
    	basisfd1  <- fd1$basis
    	basisfd2  <- fd2$basis
   		nbasis1   <- basisfd1$nbasis
    	nbasis2   <- basisfd2$nbasis
    	rangeval1 <- basisfd1$rangeval
    	rangeval2 <- basisfd2$rangeval
    	if (any(rangeval1 != rangeval2)) stop(
      		"The ranges of the arguments are not equal.")
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
    	fdnames[[3]] <- paste(fdnames1[[3]],"*",fdnames2[[3]])
    	fdquot   <- fd(coefquot, basisfd, fdnames)
    	return(fdquot)
  	} else {
	   #  first object is functional and second object is numeric
    	if (!is.numeric(fd2)) stop(
      		"second argument for * is not numeric.")
    	coef     <- fd1$coefs
    	coefd    <- dim(coef)
    	basisfd  <- fd1$basis
    	nbasis   <- basisfd$nbasis
    	rangeval <- basisfd$rangeval
    	neval    <- max(10*nbasis + 1,101)
    	evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
    	fdarray  <- eval.fd(evalarg, fd1)/fd2
    	coefquot <- project.basis(fdarray, evalarg, basisfd)
    	fdnames  <- fd1$fdnames
    	fdnames[[3]] <- paste(as.character(fd2),"*",fdnames[[3]])
    	fdquot   <- fd(coefquot, basisfd, fdnames)
    	return(fdquot)
 	}
}

#  power method for "fd"

"^.fd" <- function(fdobj, power)
{
   	#  Arguments:
	#  FD     ...  A functional data object
	#  POWER  ...  An exponent
	#  Returns:
	#  FDAPOWR  ...  A functional data object that is FD to the power POWER
	#  Last modified:  17 September 2005
	
  	if ((!(inherits(fdobj, "fd")))) stop(
      "First argument for ^ is not a functional data object.")
  	if ((!(is.numeric(power)))) stop(
      "Second argument for ^ is not numeric.")
  	coef     <- fdobj$coefs
  	coefd    <- dim(coef)
  	basisfd  <- fdobj$basis
  	nbasis   <- basisfd$nbasis
  	rangeval <- basisfd$rangeval
  	neval    <- max(10*nbasis + 1,101)
  	evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
  	fdnames  <- fdobj$fdnames
  	fdarray  <- eval.fd(evalarg, fdobj)^power
  	coefpowr <- project.basis(fdarray, evalarg, basisfd)
  	fdpowr   <- fd(coefpowr, basisfd, fdnames)
  	fdpowr
}


#  sqrt method for "fd"

sqrt.fd <- function(fdobj)
{
   	#  Arguments:
	#  FDPBJ ...  A functional data object
	#  Returns:
	#  FDAROOT  ...  A functional data object that is the square root of FDOBJ
	#  Last modified:  17 September 2005
	
  	if ((!(inherits(fdobj, "fd")))) stop(
      "First argument for ^ is not a functional data object.")
  	coef     <- fdobj$coefs
  	coefd    <- dim(coef)
  	basisfd  <- fdobj$basis
  	nbasis   <- basisfd$nbasis
  	rangeval <- basisfd$rangeval
  	neval    <- max(10*nbasis + 1,101)
  	evalarg  <- seq(rangeval[1],rangeval[2], length=neval)
  	fdnames  <- fdobj$fdnames
  	fdarray  <- sqrt(eval.fd(evalarg, fdobj))
  	coefroot <- project.basis(fdarray, evalarg, basisfd)
  	fdroot   <- fd(coefroot, basisfd, fdnames)
  	return(fdroot)
}


