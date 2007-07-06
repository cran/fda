pda.fd  <-  function(xfdlist, bwtlist=NULL, awtlist=NULL, ufdlist=NULL,
                     nfine=501)
{
#PDA computes the basis function expansions of the
#  estimates of the coefficient functions a_k(t) and b_j(t)
#  in the possibly nonhomogeneous linear differential operator
#
#    Lx(t) <- a_1(t)u_1(t) +  + a_k(t)u_K(t) +
#       b_0(t)x(t) + b_1(t)Dx(t) +  + b_{m-1}D^{m-1}x(t) + D^m x(t)
#
#  of order m <- DIFEORDER that minimizes in a least squares sense the residual
#  functions f(t) <- Lx(t).
#
#  If (DIFEORDER <- 0, PDALIST fits the varying coefficient or pointwise
#  linear model using the functions x(t) as dependent variables and
#  the forcing functions u(t) as indep}ent variables.  In this case,
#  there must be at least one forcing function.
#
#  The functions x(t) are in functional data object XFDOBJ.
#  The forcing functions u_k(t) are in functional data object UFDOBJ.
#  The coefficient functions for u_k(t) and x(t) are expanded in terms of the
#  basis functions specified in AWTLIST and BWTLIST, respectively.
#
#  The functions u_k(t) and x(t) are assumed to be vector valued
#    of dimension J.
#  That is, the differential equation can be a system of J equations rather
#    than a single equation.
#  Each coefficient function b_j(t) is matrix valued, with a column
#    for each scalar function in the system.

#  Arguments:
#  XFDLIST     list array of functional data objects for the functions
#                 whose derivatives define the DIFE
#                 dimensions are J and 1
#  BWTLIST     list array of functional parameter objects defining
#                 weight functions for functions x, Dx, ...
#                 dimensions are J, J and DIFEORDER
#  AWTLIST     list array of functional parameter objects defining
#                 weight function specs for u-variables
#                 dimensions are J and K
#  UFDLIST     list array of functional data objects for the
#                 independent variables or u-variables
#                 dimensions are J and K
#  NFINE       number of sampling points for numerical integration

#  The value in each list of XFDLIST, UFDLIST, is a
#      scalar FD object.
#  The value in each list of AWTLIST AND BWTLIST, is a
#      scalar FDPAR object.

#  The value in each list of AWTLIST and BWTLIST is a list containing:
#       FD object containing to be used as the fixed value if (not estimated
#       ESTIMATE: TRUE if (weight is to be estimated, FALSE if (not.
#  We are not bothering with smoothing at this point.

#  Returns:
#  BWTLIST     list array of weights for x functions
#                 dimensions are J, J and DIFEORDER
#  RESFDLIST   FD object for residual functions.
#  AWTLIST     list array of weights for u-variables
#                 dimension J and K

#  last modified 2007.11.28 by Spencer Graves
#  previously modified 4 December 2005

#  check dimensions of the lists

# check XFDLIST

if (inherits(xfdlist, "fd")) xfdlist = list(xfdlist)

if (!inherits(xfdlist, "list")) stop(
		"XFDLIST is neither a list or a FD object")

if (length(dim(xfdlist)) > 0) stop("XFDLIST is not a vector.")
	
nvar <- length(xfdlist)

#  ----------------------------------------------------------------
#     For efficiency, there are two versions of this code:
#     one for a single variable, and another for multiple variables.
#  ----------------------------------------------------------------

if (nvar == 1) {

#  ----------------------------------------------------------------
#                   Single variable case
#  ----------------------------------------------------------------

difeorder <- length(bwtlist)
difeordp1 <- difeorder + 1

xfdobj <- xfdlist[[1]]
xbasis <- xfdobj$basis
xcoef  <- xfdobj$coefs
xrange <- xbasis$rangeval

#  check the dimensions of UFDLIST and AWTLIST

if (is.null(ufdlist) | is.null(awtlist)) {
	nforce  <- 0
} else {
	nforce <- length(ufdlist)
	ufd1   <- ufdlist[[1]]
	urange <- ufd1$basis$rangeval
    if (length(awtlist) != nforce)
        stop("The length of AWTLIST is incorrect.")
}

#  check to see if there is anything to estimate

if (difeorder == 0 && nforce == 0)
    stop("There are no coefficient functions to estimate.")

ncurve     <- dim(xcoef)[2]

nbasmax <- 0  #  This will be the maximum number of basis functions

#  check UFDLIST and AWTLIST

if (nforce > 0) {
    for (iu in 1:nforce) {
       if (!inherits(ufdlist[[iu]], "fd")) stop(
			paste("UFDLIST[[",iu,"]] is not a functional data object."))
		ufdi     <- ufdlist[[iu]]
       urange   <- ufdi$basis$rangeval
       #  check that urange is equal to xrange
		if (any(urange != xrange)) stop(
			"XRANGE and URANGE are not identical")
		afdPari <- awtlist[[iu]]
    	afdi   <- afdPari$fd
    	if (!inherits(afdi, "fd")) stop(
			"AFDI is not a functional data object.")
    	basisi <- afdi$basis
    	if (any(basisi$rangeval != urange)) stop(
			"Ranges are incompatible for AWTLIST.")
    	nbasmax <- max(c(nbasmax,basisi$nbasis))
    }
}

#  check BWTLIST

if (difeorder > 0) {
    for (j in 1:difeorder) {
		bfdParj <- bwtlist[[j]]
       bfdj <- bfdParj$fd
       if (!inherits(bfdj, "fd")) stop(
				paste("BWTLIST[[",iu,
                    "]] is not a functional data object."))
       basisj <- bfdj$basis
       if (any(basisj$rangeval != xrange)) stop(
			"Ranges are incompatible for BWTLIST.")
       nbasmax <- max(c(nbasmax,basisj$nbasis))
    }
}

#  Set up sampling values to be used in numerical integration
#    and set up matrix of basis values.  The number of sampling
#  NFINE is here set to a usually workable value if too small.

if (nfine < 5*nbasmax) nfine <- 5*nbasmax

deltax <- (xrange[2]-xrange[1])/(nfine-1)
tx     <- seq(xrange[1],xrange[2],deltax)

if (nforce > 0) {
    deltau <- (urange[2]-urange[1])/(nfine-1)
    tu     <- seq(urange[1],urange[2],deltau)
}

#  set up  YARRAY to hold values of x functions and their derivatives

yarray <- array(0,c(nfine,ncurve,difeordp1))
for (j in 1:difeordp1) yarray[,,j] <- eval.fd(tx, xfdobj, j-1)

#  set up  UARRAY to hold values of u functions

if (nforce > 0) {
    uarray <- array(0,c(nfine,ncurve,nforce))
    for (iu in 1:nforce)
        uarray[,,iu] <- eval.fd(tu, ufdlist[[iu]])
}

#  set up array YPROD to hold mean of products of values in YARRAY

mmat  <- m2ij(nvar,difeordp1)
yprod <- array(0,c(nfine,difeordp1,difeordp1))
for (j1 in 1:difeordp1) for (j2 in 1:j1) {
        if (ncurve == 1)
			yprodval <- yarray[,1,j1]*yarray[,1,j2]
        else             yprodval <-
			apply(yarray[,,j1]*yarray[,,j2],1,mean)
        yprod[,j1,j2] <- yprodval
        yprod[,j2,j1] <- yprodval
}

#  set up array YUPROD to hold mean of u-variables u times
#    x functions and their derivatives

if (nforce > 0) {
    yuprod <- array(0,c(nfine, nforce, difeordp1))
    for (iu in 1:nforce) 
		for (j1 in 1:difeordp1) {
            if (ncurve == 1) {
                yuprodval <- yarray[,1,j1]*uarray[,iu]
            } else {
                yuprodval <-
                    apply(yarray[,,j1]*uarray[,,iu],1,mean)
            }
            yuprod[,iu,j1] <- yuprodval
    }
}

#  set up array UPROD to hold mean of products of u-variables u

if (nforce > 0) {
    uprod <- array(0,c(nfine, nforce, nforce))
    for (iu in 1:nforce) for (ju in 1:iu) {
	        if (ncurve == 1) uprodval <- uarray[,iu]*uarray[,ju]
	        else             uprodval <- apply(uarray[,,iu]*uarray[,,ju],1,mean)
            uprod[,iu,ju] <- uprodval
            uprod[,ju,iu] <- uprodval
    }
}

#  set up an index array and some arrays of 1's

onesn <- rep(1,nfine)

#  set up array to hold coefficients for basis expansions

if (nforce > 0) aarray <- matrix(0,nfine,nforce)
else            aarray <- NULL

if (difeorder > 0) barray <- matrix(0,nfine,difeorder)
else            barray <- NULL


#  --------------  beginning of loop through variables  -------------------

#  get number of coefficients to be estimated for this equation

# loop through u-variables

neqns  <- 0

if (nforce > 0) {
	for (iu in 1:nforce) {
    	afdPari <- awtlist[[iu]]
    	if (afdPari$estimate)
        	neqns <- neqns + afdPari$fd$basis$nbasis
	}
}

# loop through x functions and their derivatives

for (j1 in 1:difeorder) {
    bfdParj <- bwtlist[[j1]]
    if (bfdParj$estimate)
        neqns <- neqns + bfdParj$fd$basis$nbasis
}

if (neqns < 1) stop(
		"Number of equations to solve is not positive.")
		
#  set up coefficient array and right side array for linear equation

cmat   <- matrix(0,neqns, neqns)
dmat   <- matrix(0,neqns, 1)

#  evaluate default weight functions for this variable

if (nforce > 0) {
	for (iu in 1:nforce) {
    	afdPari     <- awtlist[[iu]]
    	aarray[,iu] <- eval.fd(tu, afdPari$fd)
	}
}

for (j1 in 1:difeorder) {
    bfdParj <- bwtlist[[j1]]
    barray[,j1] <- eval.fd(tx, bfdParj$fd)
}

#  loop through equations,
#    corresponding to rows for CMAT and DMAT

#  loop through equations for u-variables

mi12 <- 0
if (nforce > 0) {
for (iu1 in 1:nforce) {
    afdPari1   <- awtlist[[iu1]]
    if (afdPari1$estimate) {
        abasisi1    <- afdPari1$fd$basis
        abasismati1 <- getbasismatrix(tu, abasisi1)
        mi11 <- mi12 + 1
        mi12 <- mi12 + abasisi1$nbasis
        indexi1 <- mi11:mi12
        #  DMAT entry for u-variable
        weighti1 <- yuprod[,iu1,difeordp1]
        dmat[indexi1] <-
            trapzmat(abasismati1,onesn,deltax,weighti1)
        # add terms corresponding to x-derivate weights
        # that are not estimated
        for (j1 in 1:difeorder) {
	         bfdParij <- bwtlist[[j1]]
	         if (!bfdParij$estimate) {
		         weightij <- yuprod[,iu1,j1]
		         dmat[indexi1] <- dmat[indexi1] +
		             trapzmat(abasismati1, barray[,j1],
		                 deltax, weightij)
	         }
        }
        #  loop through weight functions to be estimated,
        #    corresponding to columns for CMAT
        #  begin with u-variables
        mi22 <- 0
        for (iu2 in 1:nforce) {
            afdPari2   <- awtlist[[iu2]]
            if (afdPari2$estimate) {
                abasisi2    <- afdPari2$fd$basis
                abasismati2 <- getbasismatrix(tu, abasisi2)
                weighti2    <- uprod[,iu1,iu2]
                Cprod  <- trapzmat(abasismati1, abasismati2,
                                       deltau, weighti2)
                mi21 <- mi22 + 1
                mi22 <- mi22 + abasisi2$nbasis
                indexi2 <- mi21:mi22
                #  coefficient matrix CMAT entry
                cmat[indexi1,indexi2] <- Cprod
            }
        }
        #  remaining columns:
        #    loop through u-variable -- x-derivative pairs
        mij22 <- mi22
        for (j2 in 1:difeorder) {
            bfdParj2     <- bwtlist[[j2]]
            if (bfdParj2$estimate) {
                bbasisij2    <- bfdParj2$fd$basis
                bbasismatij2 <- getbasismatrix(tx, bbasisij2)
                weightij12   <- yuprod[,iu1,j2]
                Cprod <- trapzmat(abasismati1,bbasismatij2,
                                      deltax,weightij12)
                mij21 <- mij22 + 1
                mij22 <- mij22 + bbasisij2$nbasis
                indexij2  <- mij21:mij22
                cmat[indexi1,indexij2] <- Cprod
            }
        }
        #  add roughness penalty matrix to diagonal entries
        lambdai1 <- afdPari1$lambda
        if (lambdai1 > 0) {
	         Lfdobj <- afdPari1$Lfd
	         penmat <- lambdai1*eval.penalty(abasisi1, Lfdobj)
	         cmat[indexi1,indexi1] <- cmat[indexi1,indexi1] +
	                     penmat;
        }
    }
}
}

#  loop through equations for x-derivatives

mij12 <- mi12
for (j1 in 1:difeorder) {
    bfdParj1 <- bwtlist[[j1]]
    if (bfdParj1$estimate) {
        bbasisij1    <- bfdParj1$fd$basis
        bbasismatij1 <- getbasismatrix(tx,bbasisij1)
        mij11 <- mij12 + 1
        mij12 <- mij12 + bbasisij1$nbasis
        indexij1 <- mij11:mij12
        #  DMAT entry for u-variable -- x-derivative pair
        weightij1 <- yprod[,j1,difeordp1]
        dmat[indexij1] <-
            trapzmat(bbasismatij1,onesn,deltax,weightij1)
        #  add terms corresponding to forcing functions
        #  with unestimated coefficients
		 if (nforce > 0) {
        for (iu in 1:nforce) {
	         afdPari <- awtlist[[iu]]
	         if (!afdPari$estimate) {
		         weightijk <- yuprod[,iu,j1]
		         dmat[indexij1] <- dmat[indexij1] +
		              trapzmat(bbasisij1, aarray[,iu],
		                           deltax, weightijk)
	         }
        }
        }
        #  first columns of CMAT: u-variable entries
        mi22 <- 0
        if (nforce > 0) {
        for (iu2 in 1:nforce) {
            afdPari2 <- awtlist[[iu2]]
            if (afdPari2$estimate) {
                abasisi2    <- afdPari2$fd$basis
                abasismati2 <- getbasismatrix(tx, abasisi2)
                weighti2    <- yuprod[,iu2,j1]
                Cprod <-
                    trapzmat(bbasismatij1,abasismati2,deltax,weighti2)
                mi21 <- mi22 + 1
                mi22 <- mi22 + abasisi2$nbasis
                indexi2 <- mi21:mi22
                cmat[indexij1,indexi2] <- Cprod
            }
        }
        }
        #  remaining columns: x-derivative pairs
        mij22 <- mi22
        for (j2 in 1:difeorder) {
            bfdParj2  <- bwtlist[[j2]]
            bbasisij2 <- bfdParj2$fd$basis
            if (bfdParj2$estimate) {
                bbasismatij2 <- getbasismatrix(tx, bbasisij2)
                weightij22   <- yprod[,j1,j2]
                Cprod <- trapzmat(bbasismatij1,bbasismatij2,
                                  deltax,weightij22)
                mij21 <- mij22 + 1
                mij22 <- mij22 + bbasisij2$nbasis
                indexij2 <- mij21:mij22
                cmat[indexij1,indexij2] <- Cprod
            }
        }
        # add roughness penalty matrix to diagonal entries
        lambdaj1 <- bfdParj1$lambda
        if (lambdaj1 > 0) {
	         Lfdobj <- bfdParj1$Lfd
	         penmat <- lambdaj1*eval.penalty(bbasisij1, Lfdobj)
	         cmat[indexij1,indexij1] <- cmat[indexij1,indexij1] +
	                  penmat
        }
    }
}

#  --------------  end of loop through variables  -------------------

# solve for coefficients of basis expansions

dvec <- -symsolve(cmat,dmat)

#  set up u-function weight functions

mi2 <- 0
if (nforce > 0) {
	for (iu in 1:nforce) {
    	afdPari <- awtlist[[iu]]
    	if (afdPari$estimate) {
        	mi1 <- mi2 + 1
        	mi2 <- mi2 + afdPari$fd$basis$nbasis
        	indexi <- mi1:mi2
        	afdPari$fd$coefs <- dvec[indexi]
          awtlist[[iu]] <- afdPari
    	}
	}
}

#  set up X-function derivative weight functions

mij2 <- mi2
for (j in 1:difeorder) {
    bfdParj <- bwtlist[[j]]
    if (bfdParj$estimate) {
        mij1 <- mij2 + 1
        mij2 <- mij2 + bfdParj$fd$basis$nbasis
        indexij <- mij1:mij2
        bfdParj$fd$coefs <- dvec[indexij]
        bwtlist[[j]] <- bfdParj
    }
}

#  set up residual list RESFDLIST

#  initialize with highest order derivative for this variable
resmat  <- eval.fd(tx, xfdobj, difeorder)
#  add contributions from weighted u-functions
if (nforce > 0) {
	onesncurve <- rep(1,ncurve)
	for (iu in 1:nforce) {
		afdPari  <- awtlist[[iu]]
    	aveci    <- as.vector(eval.fd(tu, afdPari$fd))
    	umati    <- eval.fd(tu, ufdlist[[iu]])
    	aumati   <- outer(aveci,onesncurve)*umati
    	resmat   <- resmat + aumati
	}
}
#  add contributions from weighted x-function derivatives
for (j in 1:difeorder) {
	bfdParj <- bwtlist[[j]]
    bmatij <- as.vector(eval.fd(tx, bfdParj$fd))
    xmatij <- eval.fd(tx, xfdobj, j-1)
    resmat <- resmat + bmatij*xmatij
}
#  set up the functional data object
resbasis <- xbasis
resfd    <- data2fd(resmat, tx, resbasis)
resfdnames      <- xfdobj$fdnames
resfdnames[[2]] <- "Residual function"
resfdnames[[3]] <- "Residual function value"
resfd$fdnames   <- resfdnames
resfdlist       <- list(resfd)

#  ----------------------------------------------------------------
#                   End of single variable case
#  ----------------------------------------------------------------

} else {

#  ----------------------------------------------------------------
#                   Multiple variable case
#  ----------------------------------------------------------------

difeorder <- dim(bwtlist)[3]
difeordp1 <- difeorder + 1

#  check the dimensions of UFDLIST and AWTLIST

if (ufdlist == NULL || awtlist == NULL) {
    nforce  <- 0
    awtlist <- NULL
} else {
	nforce <- dim(ufdlist)[2]
    if (any(dim(ufdlist) != c(nvar,nforce)))
        stop(paste("The number of rows of UFDLIST",
                   " does not match that of XFDLIST."))
    if (any(dim(awtlist) != c(nvar,nforce)))
        stop("The dimensions of AWTLIST are incorrect.")
	awtlist <- matrix("list", nvar, difeorder)
}

#  check to see if there is anything to estimate

if (difeorder == 0 && nforce == 0)
    stop("There are no coefficient functions to estimate.")

#  check the dimensions of BWTLIST

if (difeorder == 1) {
    if (any(dim(bwtlist) != c(nvar,nvar))) stop(
		"The dimensions of BWTLIST are incorrect.")
	 bwtlist <- matrix("list", c(nvar,nvar))
} else {
    if (any(dim(bwtlist) != c(nvar,nvar,difeorder))) stop(
		"The dimensions of BWTLIST are incorrect.")
	bwtlist <- array("list", c(nvar,nvar,difeorder))
}

#  check XFDLIST and extract NCURVE and XRANGE

xfd1       <- xfdlist[[1]]
xcoef1     <- xfd1$coefs
xbasis1    <- xfd1$basis
xrange1    <- xbasis1$rangeval
ncurve     <- dim(xcoef1)[2]
resfdnames <- xfd1$fdnames
for (ivar in 1:nvar) {
	xfdi    <- xfdlist[[ivar]]
	xcoefi  <- xfdi$coefs
	xbasisi <- xfdi$basis
 	xrangei <- xbasisi$rangeval
	ncurvei <- dim(xcoefi)[2]
 	if (!inherits(xfdi, "fd"))
        stop(paste("XFDLIST[[",ivar,
                "]] is not a functional data object."))
  	if (any(xrangei != xrange1)) stop(
			"Ranges are incompatible for XFDLIST.")
    if (ncurvei != ncurve) stop(
			"Number of curves is incompatible for XFDLIST.")
}

nbasmax <- 0  #  This will be the maximum number of basis functions

#  check UFDLIST and AWTLIST
#  Note:  XRANGE and URANGE are required to be identical in this version

if (nforce > 0) {
	ufd11  <- ufdlist[[1,1]]
	urange <- ufd11$basis$rangeval
    for (ivar in 1:nvar) for (iu in 1:nforce) {
			  ufdiviu <- ufdlist[[ivar,iu]]
            if (!inherits(ufdiviu, "fd")) stop(
					paste("UFDLIST[[",ivar,",",iu,
                        "]] is not a functional data object."))

            if (any(ufdiviu$basis$rangeval != urange)) stop(
					"Ranges are incompatible for UFDLIST.")
            awtfdPari <- awtlist[[ivar,iu]]
            if (!inherits(afdPari, "fdPar")) stop(
					paste("AWTFDPAR[[",ivar,",",iu,
					      "]] is not a functional data object."))
            afdi      <- awtfdPari$fd
            basisi    <- afdi$basis
            if (any(basisi$rangeval != urange)) stop(
					"Ranges are incompatible for AWTLIST.")
            nbasmax <- max(c(nbasmax,basisi$nbasis))
    }
}

#  check BWTLIST

bwtlistd <- dim(bwtlist)

if (length(bwtlistd) != 3) stop(
	"BWTLIST does not have three dimensions.")
	
if (bwtlistd[1] != nvar | bwtlistd[2] != nvar) stop(
	"First two dimensions of BWTLIST not equal to length of XFDLIST.")

if (difeorder > 0) {
    for (ivar1 in 1:nvar) for (ivar2 in 1:nvar) for (j in 1:difeorder) {
                bfdPari1i2j <- bwtlist[[ivar1,ivar2,j]]
                if (!inherits(bfdPari1i2j, "fdPar"))
                  stop(paste("BWTLIST[[",ivar1, ",",ivar2, ",",iu,
                            "]] is not a functional parameter object."))
                basisi1i2j <- bfdPari1i2j$fd$basis
                if (any(basisi1i2j$rangeval != xrange)) stop(
						"Ranges are incompatible for BWTLIST.")
                nbasmax <- max(c(nbasmax,basisi$nbasis))

    }
}

#  At this point we assume that the ranges for XFDLIST and UFDLIST
#  are the same, but this will be changed later to allow for lags.

if (nforce > 0) if (any(xrange != urange)) stop(
			"Ranges for XFDLIST and UFDLIST are not compatible.")

#  set up sampling values to be used in numerical integration
#    and set up matrix of basis values.  The number of sampling
#  NFINE is here set to a usually workable value if too small.

if (nfine < 5*nbasmax) nfine <- 5*nbasmax

deltax <- (xrange[2]-xrange[1])/(nfine-1)
tx     <- seq(xrange[1],xrange[2],deltax)

if (nforce > 0) {
    deltau <- (urange[2]-urange[1])/(nfine-1)
    tu     <- seq(urange[1],urange[2],deltau)
}

#  set up  YARRAY to hold values of x functions and their derivatives

yarray <- array(0,c(nfine,ncurve,nvar,difeordp1))
for (ivar in 1:nvar) for (j in 1:difeordp1)
	yarray[,,ivar,j] <- eval.fd(tx, xfdlist[[ivar]], j-1)

#  set up  UARRAY to hold values of u functions

if (nforce > 0) {
    uarray <- array(0,c(nfine,nforce))
    for (iu in 1:nforce)
        uarray[,iu] <- eval.fd(tu, ufdlist[[ivar,iu]])
}

#  set up array YPROD to hold mean of products of values in YARRAY

mmat  <- m2ij(nvar,difeordp1)
yprod <- array(0,c(nfine,nvar,difeordp1,nvar,difeordp1))
for (m1 in 1:nvar*difeordp1) {
    i1 <- mmat[m1,1]
    j1 <- mmat[m1,2]
    for (m2 in 1:m1) {
        i2 <- mmat[m2,1]
        j2 <- mmat[m2,2]
        if (ncurve == 1)
            yprodval <-       yarray[,,i1,j1]*yarray[,,i2,j2]
        else
            yprodval <- apply(yarray[,,i1,j1]*yarray[,,i2,j2],2,mean)
        yprod[,i1,j1,i2,j2] <- yprodval
        yprod[,i2,j2,i1,j1] <- yprodval
    }
}

#  set up array YUPROD to hold mean of u-variables u times
#    x functions and their derivatives

if (nforce > 0) {
    yuprod <- array(0,c(nfine, nvar, nforce, difeordp1))
    for (iu in 1:nforce) {
        for (i1 in 1:nvar) {
            for (j1 in 1:difeordp1) {
                if (ncurve == 1)
                    yuprodval <- yarray[,1,i1,j1]*uarray[,iu]
                else
                    yuprodval <-
                        apply(yarray[,,i1,j1]*
                              outer(uarray[,iu],onesncurve),2,mean)
                yuprod[,i1,iu,j1] <- yuprodval
            }
        }
    }
}

#  set up array UPROD to hold mean of products of u-variables u

if (nforce > 0) {
    uprod <- array(0,c(nfine, nforce, nforce))
    for (iu in 1:nforce) for (ju in 1:iu) {
            uprodval <- uarray[,iu]*uarray[,ju]
            uprod[,iu,ju] <- uprodval
            uprod[,ju,iu] <- uprodval
    }
}

#  set up an index array and some arrays of 1"s

mmat  <- m2ij(nvar,difeorder)
onesn <- rep(1,nfine)

#  set up array to hold coefficients for basis expansions

if (nforce > 0) aarray <- matrix(0,nfine,nforce)
else            aarray <- NULL

if (difeorder > 0) barray <- array(0,c(nfine,nvar,difeorder))
else            barray <- NULL

#  --------------  beginning of loop through variables  -------------------

for (ivar in 1:nvar) {

    #  get number of coefficients to be estimated for this equation

    # loop through u-variables
    neqns  <- 0
    for (iu in 1:nforce) {
        afdPari <- awtlist[[ivar,iu]]
        if (afdPari$estimate)
            neqns <- neqns + afdPari$fd$basis$nbasis
    }
    # loop through x functions and their derivatives
    for (m2 in 1:nvar*difeorder) {
        i2 <- mmat[m2,1]
        j2 <- mmat[m2,2]
        if (difeorder == 1) bfdParij <- bwtlist[[ivar,i2]]
        else             bfdParij <- bwtlist[[ivar,i2,j2]]
        if (bfdParij$estimate)
            neqns <- neqns + bfdParij$fd$basis$nbasis
    }
    if (neqns < 1)  stop(
			"Number of equations to solve is not positive.")

    #  set up coefficient array and right side array for linear equation

    cmat   <- matrix(0,neqns, neqns)
    dmat   <- matrix(0,neqns, 1)

    #  evaluate default weight functions for this variable

    for (iu in 1:nforce) {
        afdPari <- awtlist[[ivar,iu]]
        aarray[,iu] <- eval.fd(tu, afdPari$fd)
    }
    for (i in 1:nvar) for (j in 1:difeorder) {
            if (difeorder == 1) bfdParij <- bwtlist[[ivar,i]]
            else             bfdParij <- bwtlist[[ivar,i,j]]
            barray[,i,j] <- eval.fd(tx, bfdParij$fd)
    }

    #  loop through equations,
    #    corresponding to rows for CMAT and DMAT

    #  loop through equations for u-variables

    mi12 <- 0
    for (iu1 in 1:nforce) {
        afdPari1   <- awtlist[[ivar,iu1]]
        if (afdPari1$estimate) {
            abasisi1    <- afdPari1$basis
            abasismati1 <- getbasismatrix(tu, abasisi1)
            mi11 <- mi12 + 1
            mi12 <- mi12 + abasisi1$nbasis
            indexi1 <- mi11:mi12
            #  DMAT entry for u-variable
            weighti1 <- yuprod[,ivar,iu1,difeordp1]
            dmat[indexi1] <-
                trapzmat(abasismati1,onesn,deltax,weighti1)
            #  add terms corresponding to x-derivative weights
            #  that are not estimated
            for (m in 1:(nvar*difeorder)) {
	             i <- mmat(m,1)
	             j <- mmat(m,2)
	             if (difeorder == 1) bfdParij <- bwtlist[[ivar,i]]
	             else             bfdParij <- bwtlist[[ivar,i,j]]
	             if (!bfdParij$estimate) {
		              weightij <- yuprod[,ivar,iu1,j]
		              dmat[indexi1] <- dmat[indexi1] +
		                  trapzmat(abasismati1, barray[,ivar,j],
		                           deltax, weightij)
	             }
            }
            #  loop through weight functions to be estimated,
            #    corresponding to columns for CMAT
            #  begin with u-variables
            mi22 <- 0
            for (iu2 in 1:nforce) {
                afdPari2   <- awtlist[[ivar,iu2]]
                if (afdPari2$estimate) {
                    abasisi2    <- afdPari2$fd$basis
                    abasismati2 <- getbasismatrix(tu, abasisi2)
                    weighti2    <- uprod[,iu1,iu2]
                    Cprod       <- trapzmat(abasismati1, abasismati2,
                                            deltau, weighti2)
                    mi21 <- mi22 + 1
                    mi22 <- mi22 + abasisi2$nbasis
                    indexi2 <- mi21:mi22
                    #  coefficient matrix CMAT entry
                    cmat[indexi1,indexi2] <- Cprod
                }
            }
            #  remaining columns:
            #    loop through u-variable -- x-derivative pairs
            mij22 <- mi22
            for (m2 in 1:nvar*difeorder) {
                i2 <- mmat[m2,1]
                j2 <- mmat[m2,2]
                if (difeorder == 1) bfdParij2   <- bwtlist[[ivar,i2]]
                else             bfdParij2   <- bwtlist[[ivar,i2,j2]]
                if (bfdParij2$estimate) {
                    bbasisij2    <- bfdParij2$fd$basis
                    bbasismatij2 <- getbasismatrix(tx, bbasisij2)
                    weightij12   <- yuprod[,ivar,iu1,j2]
                    Cprod        <- trapzmat(abasismati1,bbasismatij2,
											   deltax,weightij12)
                    mij21 <- mij22 + 1
                    mij22 <- mij22 + bbasisij2$nbasis
                    indexij2  <- mij21:mij22
                    cmat[indexi1,indexij2] <- Cprod
                }
            }
            #  add roughness penalty matrix to diagonal entries
            lambdai1 <- afdPari1$lambda
            if (lambdai1 > 0) {
	             Lfdobj <- afdPari1$Lfd
	             penmat <- lambdai1*eval.penalty(abasisi1,Lfdobj)
	             cmat[indexi1,indexi1] <- cmat[indexi1,indexi1] +
	                     penmat
            }
        }
    }

    #  loop through equations for x-derivatives

    mij12 <- mi12
    for (m1 in 1:nvar*difeorder) {
        i1 <- mmat[m1,1]
        j1 <- mmat[m1,2]
        if (difeorder == 1)  bfdParij1 <- bwtlist[[ivar,i1]]
        else              bfdParij1 <- bwtlist[[ivar,i1,j1]]
        if (bfdParij1$estimate) {
            bbasisij1    <- bfdParij1$fd$basis
            bbasismatij1 <- getbasismatrix(tx,bbasisij1)
            mij11 <- mij12 + 1
            mij12 <- mij12 + bbasisij1$nbasis
            indexij1 <- mij11:mij12
            #  DMAT entry for u-variable -- x-derivative pair
            weightij1 <- yprod[,i1,j1,ivar,difeordp1]
            dmat[indexij1] <- trapzmat(bbasismatij1,onesn,
                                       deltax,weightij1)
            #  add terms corresponding to forcing functions
            #  with unestimated coefficients
            for (iu in 1:nforce) {
	             afdPari <- awtlist[[ivar,iu]]
	             if (!afdPari$estimate) {
		              weightijk <- yprod[,ivar,iu,j1]
		              dmat[indexij1] <- dmat[indexij1] +
		                  trapzmat(bbasismatij1,aarray[,iu],
		                           deltax,weightijk)
		         }
	         }
            #  first columns of CMAT: u-variable entries
            mi22 <- 0
            for (iu2 in 1:nforce) {
                afdPari2  <- awtlist[[ivar,iu2]]
               if (afdPari2$estimate) {
                    abasisi2    <- afdPari2$fd$basis
                    abasismati2 <- getbasismatrix(tx, abasisi2)
                    weighti2    <- yuprod[,i1,iu2,j1]
                    Cprod <- trapzmat(bbasismatij1,abasismati2,
                                      deltax,weighti2)
                    mi21 <- mi22 + 1
                    mi22 <- mi22 + abasisi2$nbasis
                    indexi2 <- mi21:mi22
                    cmat[indexij1,indexi2] <- Cprod
                }
            }
            #  remaining columns: x-derivative pairs
            mij22 <- mi22
            for (m2 in 1:nvar*difeorder) {
                i2 <- mmat[m2,1]
                j2 <- mmat[m2,2]
                if (difeorder == 1) bfdParij2 <- bwtlist[[ivar,i2]]
                else             bfdParij2 <- bwtlist[[ivar,i2,j2]]
                bbasisij2    <- bfdParij2$fd$basis
                bbasismatij2 <- getbasismatrix(tx, bbasisij2)
                weightij22   <- yprod[,i1,j1,i2,j2]
                Cprod <-
                    trapzmat(bbasismatij1,bbasismatij2,deltax,weightij22)
                if (bfdParij2$estimate) {
                    mij21 <- mij22 + 1
                    mij22 <- mij22 + bbasisij2$nbasis
                    indexij2 <- mij21:mij22
                    cmat[indexij1,indexij2] <- Cprod
                }
            }
            #  add roughness penalty terms to diagonal entries
            lambdaij1 <- bfdParij1$lambda
            if (lambdaij1 > 0) {
	             Lfdobj <- bfdParij1$Lfd
	             penmat <- lambdaij1*eval.penalty(bbasisij1,Lfdobj)
	             cmat[indexij1,indexij1] <- cmat[indexij1,indexij1] +
	                  penmat
            }
        }
    }

    dvec <- -symsolve(cmat,dmat)

    #  set up u-function weight functions

    mi2 <- 0
    for (iu in 1:nforce) {
        afdPari <- awtlist[[ivar,iu]]
        if (afdPari$estimate) {
            mi1 <- mi2 + 1
            mi2 <- mi2 + afdPari$fd$basis$nbasis
            indexi <- mi1:mi2
            afdPari$fd$coefs <- dvec[indexi]
        }
    }

    #  set up X-function derivative weight functions

    mij2 <- mi2
    for (m1 in 1:nvar*difeorder) {
        i1 <- mmat[m1,1]
        j1 <- mmat[m1,2]
        bfdParij <- bwtlist[[ivar,i1,j1]]
        if (bfdParij$estimate) {
            mij1 <- mij2 + 1
            mij2 <- mij2 + bfdParij$fd$basis$nbasis
            indexij <- mij1:mij2
            bfdParij$fd$coefs <- dvec[indexij]
            bwtlist[[ivar,i1,j1]] <- bfdParij
        }
    }
}

#  --------------  end of loop through variables  -------------------

#  set up residual list RESFDLIST

resfdlist       <- vector("list", nvar)

for (ivar in 1:nvar) {
	xfdi     <- xfdlist[[ivar]]
	resbasis <- xfdi$basis
    #  initialize with highest order derivative for this variable
    resmat  <- eval.fd(tx, xfdi, difeorder)
    #  add contributions from weighted u-functions
    if (nforce > 0) {
		onesncurve <- rep(1,ncurve)
    	for (iu in 1:nforce) {
		 	afdPari  <- awtlist[[ivar,iu]]
        	amati    <- as.vector(eval.fd(tu, afdPari$fd))
        	umati    <- eval.fd(tu, ufdlist[[ivar,iu]])
		 	if (ncurve == 1) aumati <- amati*umati
		 	else             aumati <- outer(amati,onesncurve)*umati
        	resmat   <- resmat + aumati
    	}
    }
    #  add contributions from weighted x-function derivatives
    for (m1 in 1:nvar*difeorder) {
        i1 <- mmat[m1,1]
        j1 <- mmat[m1,2]
		 bfdParij <- bwtlist[[ivar,i1,j1]]
		 bfdij    <- bfdParij$fd
		 bvecij   <- as.vector(eval.fd(tx, bfdij))
		 if (ncurve == 1) bmatij <- bvecij
		 else             bmatij <- outer(bvecij,onesncurve)
        xmatij <- eval.fd(tx, xfdi, j1-1)
        resmat <- resmat + bmatij*xmatij
    }
    #  set up the functional data object
    resfdi            <- data2fd(resmat, tx, resbasis)
    resfdnames        <- xfdi$fdnames
    resfdnames[[2]]   <- "Residual function"
    resfdnames[[3]]   <- "Residual function value"
    resfdlist[[ivar]] <- resfdi
	
}

#  ----------------------------------------------------------------
#                   End of multiple variable case
#  ----------------------------------------------------------------

}

return(list(bwtlist=bwtlist, resfdlist=resfdlist, awtlist=awtlist))

}


#  --------------------------------------------------------------

m2ij <- function(nrow,ncol) {
#M2IJ sets up a NROW*NCOL by 2 matrix of row-col indices associated
#  with a number of matrix entries row-wise
#  Example:  m2ij(2,3) produces
#     1     1
#     1     2
#     1     3
#     2     1
#     2     2
#     2     3
	nval <- nrow*ncol
	if (nval > 0)
    	mmat <- cbind(c(outer(rep(1,ncol),(1:nrow))),
              	    c(outer((1:ncol),rep(1,nrow))))
	else mmat <- NULL
	
	return(mmat)
}

