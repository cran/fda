pda.fd <- function (fd, norder, wbasisfd=basisfd, n=5*nbasis,
                    estimate=c(FALSE,rep(TRUE,nord)),
                    lambda=rep(0,nordp1), wfd0=matrix(0,1,nordp1))
{

#  A function to compute the basis function expansion of the
#    estimate of the coefficient functions w_j(t) and forcing function f(t) 
#    in the nonhomogeneous linear differential operator
#
#    Lx(t) = f(t) +
#       w_0(t)x(t) + w_1(t)Dx(t) + ... + w_{m-1}D^{m-1}x(t) + w_m(t)D^m x(t)  
#
#    of order m = NORDER that minimizes in a least squares sense the residual
#    functions Lx(t).  The functions x(t) are in functional data object FD.
#  The coefficient functions are expanded in terms of the
#    basis functions specified in WBASISFD.

#  Arguments:
#  FD        ...  functional data object
#  NORDER    ...  order of the linear differential operator, that is, the order
#                 of the highest derivative.
#  WBASISFD  ...  basis object for the forcing function and weight functions.
#  N         ...  number of sampling points for numerical integration
#  ESTIMATE  ...  logical array of length NORDER + 1, if a value is T, the
#                 corresponding coefficient function is estimated, otherwise
#                 the target value is used.  The first value applies to the 
#                 forcing function, and if F, the forcing function is 0 and
#                 linear differential operator is homogeneous
#  LAMBDA    ...  a numerical array of length NORDER + 1 containing 
#                 penalty parameters for penalizing the departure of the
#                 estimated weight functions from those defined in WFD0.
#  WFD0      ...  A specification of a functional data object that is used for
#                 those weight functions and forcing function not estimated, 
#                 or as target functions toward which the estimated weight 
#                 functions and forcing function are smoothed. 
#                 WFD0 can either be a vector of NORDER + 1 constants, or a 
#                 functional data object with the same structure as WFD 
#                 that is returned by this function.

#  Returns:
#  WFD       ...  estimated weight functional data object.  It has NORDER + 1
#                 functions.  The first is the forcing function, and 
#                 the last NORDER of these are the weight functions w_j(t) 
#                 in the linear differential operator.


#  Last modified 3 December 2001

#  check the first argument

	if (!(inherits(fd, "fd"))) stop("Argument FD not a functional data object.")

	basisfd <- getbasis(fd)
	nbasis  <- basisfd$nbasis
	rangew  <- basisfd$rangeval

#  check second argument

	nord <- as.integer(norder)
	if (nord <= 0) stop("nord must be positive.")
	if (nord > 9)  stop("nord exceeds 9.")
	nordp1 <- nord + 1
	nordp2 <- nord + 2
	
#  check argument WFD0.  If it is a vector, convert this to a 
#    constant basis functional data object

	if (!(inherits(wfd0, "fd")) && is.numeric(wfd0)) {
	  	if (length(wfd0) != nordp1) stop(
	         "WFD0 is not a vector of NORDER + 1")
    	wbasis0 <- create.constant.basis(rangew)
    	wfd0 <- create.fd(wfd0, wbasis0)
	} else {
    	stop("WFD0 is neither a vector nor a functional data object")
	}

#  check which components are to be estimated and get that number

   if (length(estimate) != nordp1) 
		stop("ESTIMATE is not a vector of NORDER + 1")
  	estimate <- as.logical(estimate)
  	ncoef    <- sum(as.numeric(estimate))

#  get the dimensions of the data in FD

  	coef   <- getcoef(fd)
  	coefd  <- dim(coef)
  	if (is.null(coefd)) ndim <- 1 else ndim <- length(coefd)
  	if (ndim == 1) {
  		ncurve <- 1
  		nvar   <- 1
  	}
  	if (ndim == 2) {
    	ncurve <- coefd[2]
    	nvar   <- 1
  	}
  	if (ndim == 3) {
  		ncurve <- coefd[2]
  		nvar   <- coefd[3]
  	}

#  check and get the characteristics of the basis to be used 

  	typew   <- wbasisfd$type
  	nbasisw <- wbasisfd$nbasis
  	rangew  <- wbasisfd$rangeval

  	if (any(rangew != basisfd$rangeval)) stop(
    	"Weight function range not equal to range in FD")

  	if (typew == "bspline") {
    	nbreaksw <- length(wbasisfd$params)
    	norderw  <- nbasisw - nbreaksw
  	}

#  set up sampling values to be used in numerical integration
#    and set up matrix of basis values

  	delta    <- (rangew[2]-rangew[1])/(n-1)
  	x        <- seq(rangew[1],rangew[2],delta)
  	basismat <- getbasismatrix(x, wbasisfd)

#  set up array to hold values of function and derivatives

  	yarray <- array(0,c(n,ncurve,nordp2))
  	yarray[,,1] <- 1

#  set up array to hold coefficients for basis expansion

  	if (nvar == 1) {
   			pdacoef <- matrix(0,nbasisw,nordp1)
  	} else {
    		pdacoef <- array(0,c(nbasisw,nordp1,nvar))
  	}

#  --------------  beginning of loop through variables  -------------------

  	for (ivar in 1:nvar) {
	   #  fill yarray with values of functions
  		if (nvar == 1) {
    		for (j in 0:nord) yarray[,,j+2] <- eval.fd(x, fd, j)
  		} else {
  			for (j in 0:nord) yarray[,,j+2] <- eval.fd(x, fd[,ivar], j)
  		}

    	mi   <- 0
    	mij  <- 0
    	Swgt <- matrix(0,n,ncoef)
    	Rwgt <- matrix(0,n,ncoef*(ncoef+1)/2)
    	for (i in 1:nordp1) {
      		if(estimate[i]) {
        		mi <- mi + 1
        		index <- (1 + (mi-1)*nbasisw):(mi*nbasisw)
				#  information for right side of linear equation
        		Swgt[,mi] <- apply(yarray[,,i]*yarray[,,nordp2],1,mean)
				#  information for left side of coefficient matrix for linear equation
        		mj <- 0
        		for (j in 1:i) {
          			if(estimate[j]) {
            			mij <- mij + 1
            			Rwgt[,mij] <- apply(yarray[,,i]*yarray[,,j],1,mean)
          			}
        		}
      		}
    	}

		#  set up left and right sides of linear equation
			
    	result <- SRsetup(ncoef, nbasisw, Swgt, Rwgt, basismat)

    	Cmat <- result[[2]]
    	Dmat <- result[[1]]

		#  modify the left and right sides if smoothing is involved
			
    	if (any(lambda > 0)) {
      		Hmat <- getbasispenalty(wbasisfd,0)
			mi   <- 0
      		for (i in 1:nordp1) {
  				if(estimate[i]) {
					mi <- mi + 1
        			index <- (1 + (mi-1)*nbasisw):(mi*nbasisw)
        			if (lambda[i] > 0) {
          				Cmat[index,index] <- Cmat[index,index] - lambda[i]*Hmat
						if (any(getcoef(wfd0[i]) != 0))
         	 				Dmat[index,1] <- Dmat[index,1] + 
												lambda[i]*inprod(wbasisfd,wfd0[i])
        			}
				}
      		}
    	}

		#  solve the equation using Choleski decomposition
		
    	dvec <- symsolve( Cmat, -Dmat )

    	#  set up the coefficient matrix

  		dmat <- matrix(0,nbasisw,nordp1)
  		mi  <- 0
  		for (i in 1:nordp1) {
    		if(estimate[i]) {
      			mi <- mi + 1
      			index <- (1 + (mi-1)*nbasisw):(mi*nbasisw)
      			dmat[,i] <- dvec[index]
    		}
  		}
  		if (nvar == 1) pdacoef <- dmat else pdacoef[,,ivar] <- dmat

	}

#  --------------  end of loop through variables  -------------------

#  set up the functional data object WFD 

  	wfdnames <- getnames(fd)
  	names(wfdnames)[2] <- "Weight functions"
  	names(wfdnames)[3] <- "Weight value"
  	wfd <- create.fd(pdacoef, wbasisfd, wfdnames)

  	return( wfd )
}

#  ------------------------------------------------------------------------

SRsetup <- function(ncoef, nbasis, Swgt, Rwgt, basismat)
{
  #  sets up coefficient matrices for basis expansion of weight functions

  	Smat <- matrix(0, ncoef*nbasis, 1)
  	Rmat <- matrix(0, ncoef*nbasis, ncoef*nbasis)
  	n  <- nrow(Swgt)
  	m1 <- ncol(basismat)
  	m  <- 0
  	one <- rep(1, nrow(basismat))
  	for (i in 1:ncoef){
    	indexi <- (1:nbasis) + (i-1)*nbasis
    	temp     <- basismat * outer(Swgt[,i], rep(1,m1))
    	temp[1,] <- temp[1,]/2
    	temp[n,] <- temp[n,]/2
    	Smat[indexi] <- crossprod(temp, one)
    	for (j in 1:i) {
      		m <- m + 1
      		indexj <- (1:nbasis) + (j-1)*nbasis
      		temp     <- basismat * outer(Rwgt[,m],rep(1,m1))
      		temp[1,] <- temp[1,]/2
      		temp[n,] <- temp[n,]/2
      		Rmat[indexi,indexj] <- crossprod(temp, basismat)
      		if (i != j) Rmat[indexj,indexi] <- Rmat[indexi,indexj]
    	}
  	}
  	return (list(Smat, Rmat) )
}
