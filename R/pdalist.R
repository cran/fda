pdalist  <-  function(xfdlist, ufdlist=NULL, awtlist=NULL, bwtlist=NULL,
                      norder=1, nvar=1, nforce=0, nfine=101)
{
#PDA computes the basis function expansions of the
#  estimates of the coefficient functions a_k(t) and b_j(t)
#  in the possibly nonhomogeneous linear differential operator
#
#    Lx(t) <- a_1(t)u_1(t) +  + a_k(t)u_K(t) +
#       b_0(t)x(t) + b_1(t)Dx(t) +  + b_{m-1}D^{m-1}x(t) + D^m x(t)
#
#  of order m <- NORDER that minimizes in a least squares sense the residual
#  functions f(t) <- Lx(t).
#
#  If (NORDER <- 0, PDALIST fits the varying coefficient or pointwise
#  linear model using the functions x(t) as dep}ent variables and
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
#  UFDLIST     list array of independent variables or u-variables
#                 dimensions are J and K
#  AWTLIST     list array of weight function specs for u-variables
#                 dimensions are J and K
#  BWTLIST     list array of weight function specs for functions x
#                 dimensions are J, J and NORDER
#  NORDER      order of the linear differential operator, that is,
#                 the order of the highest derivative.
#  NVAR        number of variables in the linear system
#  NFORCE      number of forcing functions
#  NFINE       number of sampling points for numerical integration

#  The value in each list of XFDLIST, UFDLIST, AFDLIST and BFDLIST is a
#      scalar FD object

#  The value in each list of AWTLIST and BWTLIST is a list containing:
#       FD object containing to be used as the fixed value if (not estimated
#       ESTIMATE: T if (weight is to be estimated, F if (not.
#  We are not bothering with smoothing at this point.

#  Returns:
#  BFDLIST     list array of weights for x functions
#                 dimensions are J, J and NORDER
#  RESFDLIST   FD object for residual functions.
#  AFDLIST     list array of weights for u-variables
#                 dimension J and K

#  last modified 20 February 2003

norder <- floor(norder)
if (norder < 0) stop("NORDER is negative.")
nordp1 <- norder + 1

#  check the dimensions of UFDLIST and AWTLIST

if (is.null(ufdlist) || is.null(awtlist)) {
    afdlist <- NULL
} else {
    if (length(ufdlist) != nvar*nforce) {
        stop(paste("The number of rows of UFDLIST)",
               " does not match that of XFDLIST."))
    }
    if (length(awtlist) != nvar*nforce) {
        stop("The dimensions of AWTLIST are incorrect.")
    }
	afdlist <- vector("list", nvar*norder)
}

#  check to see if there is anything to estimate

if (norder == 0 && nforce == 0) {
    stop("There are no coefficient functions to estimate.")
}

#  check the dimensions of BWTLIST

if (norder == 0) {
    bfdlist <- NULL
}

if (norder == 1) {
    if (length(bwtlist) != nvar*nvar) {
        stop("The dimensions of BWTLIST are incorrect.")
    }
	bfdlist <- vector("list", nvar*nvar*norder)
}

if (norder > 1) {
    if (length(bwtlist) != nvar*nvar*norder) {
        stop("The dimensions of BWTLIST are incorrect.")
    }
}

#  check XFDLIST and extract NCURVE and XRANGE

for (ivar in 1:nvar) {
    if (!is.fd(xfdlist[[ivar]])) {
        stop(paste("XFDLIST",ivar,
                "is not a functional data object."))
    }
    if (ivar == 1) {	    
        xrange     <- getbasis(xfdlist[[ivar]])$rangeval
        ncurve     <- dim(getcoef(xfdlist[[ivar]]))[2]
        bfdnames   <- getnames(xfdlist[[ivar]])
        resfdnames <- bfdnames
    } else {
        if (any(getbasis(xfdlist[[ivar]])$rangeval != xrange))
            stop("Ranges are incompatible for XFDLIST.")
        if (dim(getcoef(xfdlist[[ivar]]))[2] != ncurve) 
            stop("Number of curves is incompatible for XFDLIST.")
    }
}

nbasmax <- 0  #  This will be the maximum number of basis functions

#  check UFDLIST and extract URANGE

if (nforce > 0) {
    for (ivar in 1:nvar) {
        for (iu in 1:nforce) {
            if (!is.fd(ufdlist[[(iu-1)*nforce+ivar]])) {
                stop(paste("UFDLIST{",ivar,",",iu,
                        "is not a functional data object."))
            }
            if (ivar == 1 && iu == 1) {
                urange <- getbasis(ufdlist[[(iu-1)*nforce+ivar]])$rangeval
                afdnames <- getnames(ufdlist[[(iu-1)*nforce+ivar]])
            } else {
                if (any(getbasis(ufdlist[[(iu-1)*nforce+ivar]])$rangeval != urange))
                    stop("Ranges are incompatible for UFDLIST.")
            }
        }
    }

    #  check AWTLIST and extract the max. no. basis fns.

    for (ivar in 1:nvar) {
        for (iu in 1:nforce) {
            awtlisti <- awtlist[[(iu-1)*nforce+ivar]]
            afdi       <- awtlisti[[1]]
            if (!is.fd(afdi)) {
                stop("AFDI is not a functional data object.")
            }
            basisi <- getbasis(afdi)
            if (any(basisi$rangeval != urange)) {
                stop("Ranges are incompatible for AWTLIST.")
            }
            nbasmax <- max(c(nbasmax,basisi$nbasis))
        }
    }

}

#  check BWTLIST

if (norder > 0) {
    for (ivar1 in 1:nvar) {
        for (ivar2 in 1:nvar) {
            for (j in 1:norder) {
                if (norder == 1) {
                    blist12 <- bwtlist[[(ivar2-1)*nvar+ivar1]]
                } else {
                    blist12 <-
                        bwtlist[[(j-1)*norder*nvar+(ivar2-1)*nvar+ivar1]]
                }
                if (!is.fd(blist12[[1]])) {
                    stop(paste("BWTLIST",ivar1, ",",ivar2, ",",iu,
                            "is not a functional data object."))
                }
                basisi <- getbasis(blist12[[1]])
                if (any(basisi$rangeval != xrange)) {
                    stop("Ranges are incompatible for BWTLIST.")
                }
                nbasmax <- max(c(nbasmax,basisi$nbasis))
            }
        }
    }
}

#  At this point we assume that the ranges for XFDLIST and UFDLIST
#  are the same, but this will be changed later to allow for lags.

if (nforce > 0) {
    if (any(xrange != urange)) 
        stop("Ranges for XFDLIST and UFDLIST are not compatible.")
}

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

yarray <- array(0,c(nfine,ncurve,nvar,nordp1))
for (ivar in 1:nvar) {
    for (j in 1:nordp1) {
        yarray[,,ivar,j] <- eval.fd(tx, xfdlist[[ivar]], j-1)
    }
}

#  set up  UARRAY to hold values of u functions

if (nforce > 0) {
    uarray <- array(0,c(nfine,nforce))
    for (iu in 1:nforce) {
        uarray[,iu] <- eval.fd(tu, ufdlist[[(iu-1)*nforce+ivar]])
    }
}

#  set up array YPROD to hold mean of products of values in YARRAY

mmat  <- m2ij(nvar,nordp1)
yprod <- array(0,c(nfine,nvar,nordp1,nvar,nordp1))
for (m1 in 1:nvar*nordp1) {
    i1 <- mmat[m1,1]
    j1 <- mmat[m1,2]
    for (m2 in 1:m1) {
        i2 <- mmat[m2,1]
        j2 <- mmat[m2,2]
        if (ncurve == 1) {
            yprodval <- yarray[,1,i1,j1]*yarray[,1,i2,j2]
        } else {
            yprodval <- apply(yarray[,,i1,j1]*yarray[,,i2,j2],2,mean)
        }
        yprod[,i1,j1,i2,j2] <- yprodval
        yprod[,i2,j2,i1,j1] <- yprodval
    }
}

#  set up array YUPROD to hold mean of u-variables u times
#    x functions and their derivatives

onesncurve <- rep(1,ncurve)
if (nforce > 0) {
    yuprod <- array(0,c(nfine, nvar, nforce, nordp1))
    for (iu in 1:nforce) {
        for (i1 in 1:nvar) {
            for (j1 in 1:nordp1) {
                if (ncurve == 1) {
                    yuprodval <- yarray[,1,i1,j1]*uarray[,iu]
                } else {
                    yuprodval <-
                        apply(yarray[,,i1,j1]*
                              outer(uarray[,iu],onesncurve),2,mean)
                }
                yuprod[,i1,iu,j1] <- yuprodval
            }
        }
    }
}

#  set up array UPROD to hold mean of products of u-variables u

if (nforce > 0) {
    uprod <- array(0,c(nfine, nforce, nforce))
    for (iu in 1:nforce) {
        for (ju in 1:iu) {
            uprodval <- uarray[,iu]*uarray[,ju]
            uprod[,iu,ju] <- uprodval
            uprod[,ju,iu] <- uprodval
        }
    }
}

#  set up an index array and some arrays of 1"s

mmat  <- m2ij(nvar,norder)
onesn <- rep(1,nfine)

#  set up array to hold coefficients for basis expansions

if (nforce > 0) {
    aarray <- matrix(0,nfine,nforce)
} else {
    aarray <- NULL
}

if (norder > 0) {
    barray <- array(0,c(nfine,nvar,norder))
} else {
    barray <- NULL
}


#  --------------  beginning of loop through variables  -------------------

for (ivar in 1:nvar) {

    #  get number of coefficients to be estimated for this equation

    # loop through u-variables
    neqns  <- 0
    for (iu in 1:nforce) {
        alisti <- awtlist[[(iu-1)*nforce+ivar]]
        if (alisti[[2]]) {
            neqns <- neqns + getbasis(alisti[[1]])$nbasis
        }
    }
    # loop through x functions and their derivatives
    for (m2 in 1:nvar*norder) {
        i2 <- mmat[m2,1]
        j2 <- mmat[m2,2]
        if (norder == 1) {
            blistij <- bwtlist[[(i2-1)*nvar+ivar]]
        } else {
            blistij <- bwtlist[[(j2-1)*norder*nvar+(i2-1)*nvar+ivar]]
        }
        if (blistij[[2]]) {
            neqns <- neqns + getbasis(blistij[[1]])$nbasis
        }
    }
    if (neqns < 1) stop("Number of equations to solve is not positive.")

    #  set up coefficient array and right side array for linear equation

    cmat   <- matrix(0,neqns, neqns)
    dmat   <- matrix(0,neqns, 1)

    #  evaluate default weight functions for this variable

    for (iu in 1:nforce) {
        alisti <- awtlist[[(iu-1)*nforce+ivar]]
        aarray[,iu] <- eval.fd(tu, alisti[[1]])
    }
    for (i in 1:nvar) {
        for (j in 1:norder) {
            if (norder == 1) {
                blistij <- bwtlist[[(i-1)*nvar+ivar]]
            } else {
                blistij <- bwtlist[[(j-1)*norder*nvar+(i-1)*nvar+ivar]]
            }
            barray[,i,j] <- eval.fd(tx,blistij[[1]])
        }
    }

    #  loop through equations,
    #    corresponding to rows for CMAT and DMAT

    #  loop through equations for u-variables

    mi12 <- 0
    for (iu1 in 1:nforce) {
        alisti1   <- awtlist[[(iu1-1)*nvar+ivar]]
        if (alisti1[[2]]) {
            abasisi1    <- getbasis(alisti1[[1]])
            abasismati1 <- getbasismatrix(tu, abasisi1)
            mi11 <- mi12 + 1
            mi12 <- mi12 + abasisi1$nbasis
            indexi1 <- mi11:mi12
            #  DMAT entry for u-variable
            weighti1 <- yuprod[,ivar,iu1,nordp1]
            dmat[indexi1] <-
                trapzmat(abasismati1,onesn,deltax,weighti1)
            #  add terms corresponding to x-derivative weights
            #  that are not estimated
            for (m in 1:nvar*norder) {
                i <- mmat[m,1]
                j <- mmat[m,2]
                blistij <- bwtlist[[(j-1)*nvar*norder+(i-1)*nvar+ivar]]
                if (!blistij[[2]]) {
                    weightij <- yuprod[,ivar,k1,j]
                    dmat(indexi1) <- dmat(indexi1) + 
                        trapzmat(abasismatk1,barray[,ivar,j],deltax,weightij)
                }
            }
            #  loop through weight functions to be estimated,
            #    corresponding to columns for CMAT
            #  begin with u-variables
            mi22 <- 0
            for (iu2 in 1:nforce) {
                alisti2   <- awtlist[[(iu2-1)*nvar+ivar]]
                if (alisti2[[2]]) {
                    abasisi2    <- getbasis(alisti2[[1]])
                    abasismati2 <- getbasismatrix(tu, abasisi2)
                    weighti2    <- uprod[,iu1,iu2]
                    Cprod  <-
                        trapzmat(abasismati1, abasismati2, deltau, weighti2)
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
            for (m2 in 1:nvar*norder) {
                i2 <- mmat[m2,1]
                j2 <- mmat[m2,2]
                blistij2     <- bwtlist[[(j2-1)*norder*nvar+(i2-1)*nvar+ivar]]
                bbasisij2    <- getbasis(blistij2[[1]])
                bbasismatij2 <- getbasismatrix(tx, bbasisij2)
                weightij12   <- yuprod[,ivar,iu1,j2]
                Cprod <-
                    trapzmat(abasismati1,bbasismatij2,deltax,weightij12)
                mij21 <- mij22 + 1
                mij22 <- mij22 + bbasisij2$nbasis
                indexij2  <- mij21:mij22
                cmat[indexi1,indexij2] <- Cprod
            }
            #  add roughness penalty matrix to diaginal entries
            if (alisti1[[3]] > 0.0) {
                Lfdobj <- alistk2[[4]]
                penmat <- alistk2[[3]]*eval.penalty(abasisk1, Lfdobj)
                cmat[indexi1,indexi1] <- cmat[indexi1,indexi1] + penmat
            }
        }
    }

    #  loop through equations for x-derivatives

    mij12 <- mi12
    for (m1 in 1:nvar*norder) {
        i1 <- mmat[m1,1]
        j1 <- mmat[m1,2]
        blistij1 <- bwtlist[[(j1-1)*norder*nvar+(i1-1)*nvar+ivar]]
        if (blistij1[[2]]) {
            bbasisij1    <- getbasis(blistij1[[1]])
            bbasismatij1 <- getbasismatrix(tx,bbasisij1)
            mij11 <- mij12 + 1
            mij12 <- mij12 + bbasisij1$nbasis
            indexij1 <- mij11:mij12
            #  DMAT entry for u-variable -- x-derivative pair
            weightij1 <- yprod[,i1,j1,ivar,nordp1]
            dmat[indexij1] <-
                trapzmat(bbasismatij1,onesn,deltax,weightij1)
            #  add terms corresponding to forcing functions with
            #  unestimated coefficients
            for (k in 1:nforce) {
                alistk <- awtlist[[(k-1)*nvar+ivar]]
                if (!alistk[[2]]) {
                    weightijk      <- yuprod[,ivar,k,j1]
                    dmat[indexij1] <- dmat[indexij1] + 
                        trapzmat(bbasismatij1,aarray[,k],deltax,weightijk)
                }
            }
            #  first columns of CMAT: u-variable entries
            mi22 <- 0
            for (iu2 in 1:nforce) {
                alisti2 <- awtlist[[(iu2-1)*nvar+ivar]]
                if (alisti2[[2]]) {
                    abasisi2    <- getbasis(alisti2[[1]])
                    abasismati2 <- getbasismatrix(tx, abasisi2)
                    weighti2    <- yuprod[,i1,iu2,j1]
                    Cprod <-
                        trapzmat(bbasismatij1,abasismati2,deltax,weighti2)
                    mi21 <- mi22 + 1
                    mi22 <- mi22 + abasisi2$nbasis
                    indexi2 <- mi21:mi22
                    cmat[indexij1,indexi2] <- Cprod
                }
            }
            #  remaining columns: x-derivative pairs
            mij22 <- mi22
            for (m2 in 1:nvar*norder) {
                i2 <- mmat[m2,1]
                j2 <- mmat[m2,2]
                blistij2     <- bwtlist[[(j2-1)*norder*nvar+(i2-1)*nvar+ivar]]
                bbasisij2    <- getbasis(blistij2[[1]])
                bbasismatij2 <- getbasismatrix(tx, bbasisij2)
                weightij22   <- yprod[,i1,j1,i2,j2]
                Cprod <-
                    trapzmat(bbasismatij1,bbasismatij2,deltax,weightij22)
                mij21 <- mij22 + 1
                mij22 <- mij22 + bbasisij2$nbasis
                indexij2 <- mij21:mij22
                cmat[indexij1,indexij2] <- Cprod
            }
            #  add roughness penalty matrix to diagonal entries
            if (blistij1[[3]] > 0.0) {
                Lfdobj <- blistij1[[4]]
                penmat <- blistij1[[3]]*eval.penalty(bbasisij1, Lfdobj)
                cmat[indexij1,indexij1] <- cmat[indexij1,indexij1] + penmat
            }
        }
    }


    dvec <- -solve(cmat,dmat)

    #  set up u-function weight functions

    mi2 <- 0
    for (iu in 1:nforce) {
        alisti <- awtlist[[(iu-1)*nforce+ivar]]
        if (alisti[[2]]) {
            mi1 <- mi2 + 1
            mi2 <- mi2 + getbasis(alisti[[1]])$nbasis
            indexi <- mi1:mi2
            afdlist[[(iu-1)*nforce+ivar]] <- putcoef(dvec[indexi], alisti[[1]])
        } else {
            afdlist[[(iu-1)*nforce+ivar]] <- alisti[[1]]
        }
    }

    #  set up X-function derivative weight functions

    mij2 <- mi2
    for (m1 in 1:nvar*norder) {
        i1 <- mmat[m1,1]
        j1 <- mmat[m1,2]
        blistij <- bwtlist[[(j1-1)*norder*nvar+(i1-1)*nvar+j1]]
        if (blistij[[2]]) {
            mij1 <- mij2 + 1
            mij2 <- mij2 + getbasis(blistij[[1]])$nbasis
            indexij <- mij1:mij2
            bfdlist[[(j1-1)*nvar*nvar+(i1-1)*nvar+j1]] <- 
                                  putcoef(dvec[indexij], blistij[[1]])
        } else {
            bfdlist[[(j1-1)*nvar*nvar+(i1-1)*nvar+j1]] <- blistij[[1]]
        }
    }

}

#  --------------  end of loop through variables  -------------------

#  set up residual list RESFDLIST

resfdlist       <- list(nvar)
resfdnames[[2]] <- "Residual function"
resfdnames[[3]] <- "Residual function value"

resbasis <- getbasis(xfdlist[[1]])
for (ivar in 1:nvar) {
    #  initialize with highest order derivative for this variable
    resmat  <- eval.fd(tx, xfdlist[[ivar]], norder)
    #  add contributions from weighted u-functions
    for (iu in 1:nforce) {
        amati    <- eval.fd(tu, afdlist[[(iu-1)*nforce+ivar]])
        umati    <- eval.fd(tu, ufdlist[[(iu-1)*nforce+ivar]])
		 if (ncurve == 1) {
			aumati <- amati*umati
		 } else {
        	aumati   <- outer((amati*umati),onesncurve)
		 }
        resmat   <- resmat + aumati
    }
    #  add contributions from weighted x-function derivatives
    for (m1 in 1:nvar*norder) {
        i1 <- mmat[m1,1]
        j1 <- mmat[m1,2]
		 if (ncurve == 1) {
			bmatij <- eval.fd(tx, bfdlist[[(j1-1)*nvar*nvar+(i1-1)*nvar+j1]])
		 } else {
        	bmatij <- outer(
                    eval.fd(tx, bfdlist[[(j1-1)*nvar*nvar+(i1-1)*nvar+j1]]),
                    onesncurve)
        }
        xmatij <- eval.fd(tx, xfdlist[[i1]], j1-1)
        resmat <- resmat + bmatij*xmatij
    }
    #  set up the functional data object
    resfdi <- data2fd(resmat, tx, resbasis)
    resfdi$names <- resfdnames
    resfdlist[[ivar]] <- resfdi
}

invisible(list(bfdlist, resfdlist, afdlist))

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
	if (nval > 0) {
    	mmat <- cbind(c(outer(rep(1,ncol),(1:nrow))),
              	    c(outer((1:ncol),rep(1,nrow))))
	} else {
    	mmat <- NULL
	}
	return(mmat)
}
