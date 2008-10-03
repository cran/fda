fRegress <- function(yfdPar, xfdlist, betalist, wt=rep(1,N)) {

#  FREGRESS  Fits a functional linear model using multiple
#  functional independent variables with the dependency being
#  pointwise or concurrent.
#  The case of a scalar independent variable is included by treating
#  it as a functional independent variable with a constant basis
#  and a unit coefficient.
#
#  Arguments:
#  YFDPAR   ... an object for the dependent variable,
#               which may be:
#                   a functional data object,
#                   a functional parameter (fdPar) object, or
#                   a vector
#  XFDLIST  ... a list object of length p with each list
#               containing an object for an independent variable.
#               the object may be:
#                   a functional data object or
#                   a vector
#               if XFDLIST is a functional data object or a vector,
#               it is converted to a list of length 1.
#  BETALIST ... a list object of length p with each list
#               containing a functional parameter object for
#               the corresponding regression function.  If any of
#               these objects is a functional data object, it is 
#               converted to the default functional parameter object.
#               if BETALIST is a functional parameter object
#               it is converted to a list of length 1.
#  WT       ... a vector of nonnegative weights for observations
#
#  Returns:
#  FREGRESSLIST  ... A list containing seven members with names:
#    yfdPar      ... first  argument of FREGRESS
#    xfdlist     ... second argument of FREGRESS
#    betalist    ... third  argument of FREGRESS
#    betaestlist ... estimated regression functions
#    yhatfdobj   ... functional data object containing fitted functions
#    Cmatinv     ... inverse of the coefficient matrix, needed for
#                    function FREGRESS.STDERR that computes standard errors
#    wt          ... weights for observations
#    df          ... degrees of freedom for fit

#  Last modified 23 March 2007

#  first check of XFDLIST

if (inherits(xfdlist, "fd") || inherits(xfdlist, "numeric")) 
    xfdlist <- list(xfdlist)

if (!inherits(xfdlist, "list")) stop(
	"Argument XFDLIST is not a list object.")

#  get number of independent variables p

p <- length(xfdlist)

#  check contents of XFDCELL

xerror <- FALSE
for (j in 1:p) {
	xfdj <- xfdlist[[j]]
	if (!(inherits(xfdj, "fd") || is.numeric(xfdj))) {
		print(paste("XFDLIST[[",j,"]] is neither a FD nor a NUMERIC object."))
		xerror <- TRUE
	}
}

#  first check of BETALIST

if (inherits(betalist, "fd")) betalist <- list(betalist)

if (!inherits(betalist, "list")) stop(
	"Argument BETALIST is not a list object.")

if (length(betalist) != p)  {
	cat(paste("\nNumber of regression coefficients does not match\n",
		       "number of independent variables."))
	stop("")
}

berror <- FALSE
for (j in 1:p) {
	betafdParj <- betalist[[j]]
	if (inherits(betafdParj, "fd")) {
		betafdParj    <- fdPar(betafdParj)
		betalist[[j]] <- betafdParj
	}
	if (!inherits(betafdParj, "fdPar")) {
		print(paste("BETALIST[[",j,"]] is not a FDPAR object."))
		berror <- TRUE
	}
}

if (xerror || berror) stop("")

if (inherits(yfdPar, "fdPar") || inherits(yfdPar, "fd")) {
	
    #  ----------------------------------------------------------------
    #           YFDPAR is functional for a functional parameter
    #  ----------------------------------------------------------------

    # convert to a fdaPar object if required

    if (inherits(yfdPar, "fd")) yfdPar <- fdPar(yfdPar)

    #  extract dependent variable information

    yfdobj    <- yfdPar$fd
    ylambda   <- yfdPar$lambda
    yLfdobj   <- yfdPar$Lfd
    ycoef     <- yfdobj$coefs
    ycoefdim  <- dim(ycoef)
    N         <- ycoefdim[2]
    ybasisobj <- yfdobj$basis
    rangeval  <- ybasisobj$rangeval
    ynbasis   <- ybasisobj$nbasis

    if (length(ycoefdim) > 2) stop("YFDOBJ from YFDPAR is not univariate.")

    #  check each xfdlist member.  If the object is a vector of length N,
    #  it is converted to a functional data object with a
    #  constant basis

    onebasis <- create.constant.basis(rangeval)
    onesfd   <- fd(1,onebasis)

    xerror <- FALSE
    for (j in 1:p) {
        xfdj <- xfdlist[[j]]
        if (inherits(xfdj,"fd")) {
            xcoef <- xfdj$coefs
            if (length(dim(xcoef)) > 2) stop(
					"Covariate is not univariate.")
            rangevalx  <- xfdj$basis$rangeval
            if (any(rangevalx != rangeval))  stop(
					"Range for covariate does not match that of YFDOBJ")
        } else {
            xfdlist[[j]] <- fd(matrix(xfdj,1,N), onebasis)
        }

        #  check size of coefficient array
        xfdj  <- xfdlist[[j]]
        coefj <- xfdj$coefs
        Nj <- dim(coefj)[2]
        if (Nj != N) {
            print(
               paste("Incorrect number of replications in XFDLIST",
                     "for covariate",j))
            xerror = TRUE
        }
    }
    if (xerror) stop("")

    #  check weights

    if (length(wt) != N) stop("Number of weights not equal to N.")
    if (any(wt < 0))     stop("Negative weights found.")
    if (length(wt) == 1 || var(wt) > 0) {
       wtconstant <- FALSE 
    } else {
       wtconstant <- TRUE
    }

    #  -----------------------------------------------------------
    #          set up the linear equations for the solution
    #  -----------------------------------------------------------

    #  compute the total number of coefficients to be estimated

    ncoef <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        if (betafdParj$estimate) {
        	ncoefj     <- betafdParj$fd$basis$nbasis
        	ncoef      <- ncoef + ncoefj
        }
    }

    Cmat <- matrix(0,ncoef,ncoef)
    Dmat <- rep(0,ncoef)

    #  loop through rows of CMAT

    mj2 <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        if (betafdParj$estimate) {
            betafdj    <- betafdParj$fd
            betabasisj <- betafdj$basis
            ncoefj     <- betabasisj$nbasis
            #  row indices of CMAT and DMAT to fill
            mj1    <- mj2 + 1
            mj2    <- mj2 + ncoefj
            indexj <- mj1:mj2
            #  compute right side of equation DMAT
            xfdj <- xfdlist[[j]]
            if (wtconstant) {
                xyfdj <- xfdj*yfdobj
            } else {           
                xyfdj <- (xfdj*wt)*yfdobj
            }
            wtfdj <- sum(xyfdj)
            Dmatj <- inprod(betabasisj,onesfd,0,0,rangeval,wtfdj)
            Dmat[indexj] <- Dmatj
            #  loop through columns of CMAT
            mk2 <- 0
            for (k in 1:j) {
                betafdPark <- betalist[[k]]
                if (betafdPark$estimate) {
                    betafdk    <- betafdPark$fd
                    betabasisk <- betafdk$basis
                    ncoefk     <- betabasisk$nbasis
                    #  column indices of CMAT to fill
                    mk1 <- mk2 + 1
                    mk2 <- mk2 + ncoefk
                    indexk <- mk1:mk2
                    #  set up two weight functions
                    xfdk <- xfdlist[[k]]
                    if (wtconstant) {
                        xxfdjk <- xfdj*xfdk
                    } else {           
                        xxfdjk <- (xfdj*wt)*xfdk
                    }
                    wtfdjk <- sum(xxfdjk)
                    Cmatjk <- inprod(betabasisj, betabasisk, 0, 0, rangeval, wtfdjk)
                    Cmat[indexj,indexk] <- Cmatjk
                    Cmat[indexk,indexj] <- t(Cmatjk)
                }
            }
            #  attach penalty term to diagonal block
            lambda <- betafdParj$lambda
            if (lambda > 0) {
                Lfdj  <- betafdParj$Lfd
                Rmatj <- eval.penalty(betafdParj$fd$basis, Lfdj)
                Cmat[indexj,indexj] <- Cmat[indexj,indexj] + lambda*Rmatj
            }
        }
    }

    Cmat    <- (Cmat+t(Cmat))/2

    #  check Cmat for singularity

    eigchk(Cmat)

    #  solve for coefficients defining BETA

    Lmat    <- chol(Cmat)
    Lmatinv <- solve(Lmat)
    Cmatinv <- Lmatinv %*% t(Lmatinv)

    betacoef <- Cmatinv %*% Dmat	

    #  set up fdPar objects for reg. fns. in BETAESTLIST

    betaestlist <- betalist
    mj2 <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        if (betafdParj$estimate) {
            betafdj    <- betafdParj$fd
            ncoefj     <- betafdj$basis$nbasis
            mj1    <- mj2 + 1
            mj2    <- mj2 + ncoefj
            indexj <- mj1:mj2
            coefj  <- betacoef[indexj]
            betafdj$coefs <- as.matrix(coefj)
            betafdParj$fd <- betafdj
        }
        betaestlist[[j]] <- betafdParj
    }

    #  set up fd objects for predicted values in YHATFDOBJ

    nfine     <- max(501,10*ynbasis+1)
    tfine     <- seq(rangeval[1], rangeval[2], len=nfine)
    yhatmat <- matrix(0,nfine,N)
    for (j in 1:p) {
	  xfdj       <- xfdlist[[j]]
        xmatj      <- eval.fd(tfine, xfdj)
        betafdParj <- betaestlist[[j]]
        betafdj    <- betafdParj$fd
        betavecj   <- eval.fd(tfine, betafdj)
        yhatmat    <- yhatmat + xmatj*as.vector(betavecj)
    }
    yhatfdobj <- data2fd(yhatmat, tfine, ybasisobj)

    #  set up output list object

    fRegressList <-
       	list(yfdPar      = yfdPar,
                 xfdlist     = xfdlist,
                 betalist    = betalist,
                 betaestlist = betaestlist,
                 yhatfdobj   = yhatfdobj,
                 Cmatinv     = Cmatinv,
                 wt          = wt)

 }
 else if (inherits(yfdPar,"numeric")) {

    #  ----------------------------------------------------------------
    #                   YFDPAR is scalar or multivariate
    #  ----------------------------------------------------------------

    ymat <- as.matrix(yfdPar)
    N    <- dim(ymat)[1]

    #  check weights

    if (length(wt) != N) stop("Number of weights not equal to N.")
    if (any(wt < 0))     stop("Negative weights found.")

    Zmat  <- NULL
    Rmat  <- NULL
    pjvec <- rep(0,p)
    ncoef <- 0
    for (j in 1:p) {
        xfdj       <- xfdlist[[j]]
        if (inherits(xfdj, "fd")) {
            xcoef  <- xfdj$coefs
            Nj 	 <- dim(xcoef)[2]
            if (Nj != N) stop(
				"Coefficient matrix has the wrong number of columns.")
            xbasis     <- xfdj$basis
            betafdParj <- betalist[[j]]
            bbasis     <- betafdParj$fd$basis
            bnbasis    <- bbasis$nbasis
            pjvec[j]   <- bnbasis
            Jpsithetaj <- inprod(xbasis,bbasis)
            Zmat       <- cbind(Zmat,crossprod(xcoef,Jpsithetaj)) 
            if (betafdParj$estimate) {
                lambdaj    <- betafdParj$lambda
                if (lambdaj > 0) {
                    Lfdj   <- betafdParj$Lfd
				       Rmatj <- lambdaj*eval.penalty(bbasis,Lfdj)
                }
                else Rmatj <- matrix(0,bnbasis,bnbasis)
                if (ncoef > 0) {
	                zeromat <- matrix(0,ncoef,bnbasis)
                    Rmat  <- rbind(cbind(Rmat,       zeromat),
                                   cbind(t(zeromat), Rmatj))
                } else Rmat <- Rmatj
                ncoef <- ncoef + bnbasis
            }
        }
        else if (inherits(xfdj, "double")) {
            Zmatj  <- xfdj
            Nj     <- dim(Zmatj)[1]
            if (Nj != N) stop(
				paste("Matrix in XFDLIST[[",j,"has the wrong number of rows."))
            Zmat   <- cbind(Zmat,Zmatj)
            ncoefj <- dim(Zmatj)[2]
            Rmatj  <- matrix(0,ncoefj,ncoefj)
            Rmat   <- rbind(cbind(Rmat, matrix(0,ncoef,ncoefj)),
                           cbind(matrix(0,ncoefj,ncoef), Rmatj))
            ncoef  <- ncoef + ncoefj
          }
        else stop(
			paste("XFDLIST[[",j,
			      "is neither a functional nor a multivariate object."))

    }

    #  -----------------------------------------------------------
    #          set up the linear equations for the solution
    #  -----------------------------------------------------------

    #  solve for coefficients defining BETA

    if (any(wt != 1)) {
        rtwt   <- sqrt(wt)
        Zmatwt <- Zmat*rtwt
        ymatwt <- ymat*rtwt
        Cmat   <- t(Zmatwt) %*% Zmatwt + Rmat
        Dmat   <- t(Zmatwt) %*% ymatwt
    } else {
        Cmat <- t(Zmat) %*% Zmat + Rmat
        Dmat <- t(Zmat) %*% ymat
    }

    eigchk(Cmat)

    Cmatinv  <- solve(Cmat)
	 
    betacoef <- Cmatinv %*% Dmat


    #  compute and print degrees of freedom measure

    df <- sum(diag(Zmat %*% Cmatinv %*% t(Zmat)))

    #  set up fdPar object for BETAESTFDPAR

    betaestlist <- betalist
    onebasis    <- create.constant.basis(c(0,1))
    mj2 <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        betafdj    <- betafdParj$fd
        ncoefj     <- betafdj$basis$nbasis
        mj1    <- mj2 + 1
        mj2    <- mj2 + ncoefj
        indexj <- mj1:mj2
        betacoefj <- betacoef[indexj]
        if (inherits(xfdj, "fd")) {
            betaestfdj       <- betafdj
            betaestfdj$coefs <- as.matrix(betacoefj)
            betaestfdParj    <- betafdParj
            betaestfdParj$fd <- betaestfdj
            betaestlist[[j]] <- betaestfdParj
        }
        else{
            betaestfdj       <- fd(t(betacoefj),onebasis)
            betaestfdParj    <- betafdParj
            betaestfdParj$fd <- betaestfdj
            betaestlist[[j]] <- betaestfdParj
        }
    }

    #  set up fd object for predicted values

    yhatmat <- matrix(0,N,1)
    for (j in 1:p) {
        xfdj <- xfdlist[[j]]
        if (inherits(xfdj, "fd")) {
            xbasis     <- xfdj$basis
            xnbasis    <- xbasis$nbasis
            xrng       <- xbasis$rangeval
            nfine      <- max(501,10*xnbasis+1)
            tfine      <- seq(xrng[1], xrng[2], len=nfine)
            deltat     <- tfine[2]-tfine[1]
            xmat       <- eval.fd(tfine, xfdj)
            betafdParj <- betaestlist[[j]]
            betafdj    <- betafdParj$fd
            betamat    <- eval.fd(tfine, betafdj)
            fitj       <- deltat*
				(crossprod(xmat,betamat) -
						0.5*(outer(xmat[1,    ],betamat[1,    ]) +
				            outer(xmat[nfine,],betamat[nfine,])))
            yhatmat    <- yhatmat + fitj
        }
        else{
	       betaestfdParj <- betaestlist[[j]]
            betavecj      <- betaestfdParj$fd$coefs
            yhatmat       <- yhatmat + xfdj %*% t(betavecj)
        }
    }
    yhatfdobj <- yhatmat

    #  set up output list object

    fRegressList <-
       	list(yfdPar      = yfdPar,
                 xfdlist     = xfdlist,
                 betalist    = betalist,
                 betaestlist = betaestlist,
                 yhatfdobj   = yhatfdobj,
                 Cmatinv     = Cmatinv,
                 wt          = wt,
                 df          = df)
}
 else {
    	#  YFDOBJ is neither functional nor multivariate
    	stop("YFDOBJ is neither functional nor multivariate.")
	fRegressList <- list(0)
}


return(fRegressList)

}

#  ------------------------------------------------------------------------

eigchk <- function(Cmat) {
	
    #  check Cmat for singularity

    eigval <- eigen(Cmat)$values
    ncoef  <- length(eigval)
    if (eigval[ncoef] < 0) {
	     neig <- min(length(eigval),10)
        cat("\nSmallest eigenvalues:\n")
        print(eigval[(ncoef-neig+1):ncoef])
        cat("\nLargest  eigenvalues:\n")
        print(eigval[1:neig])
        stop("Negative eigenvalue of coefficient matrix.")
    }
    if (eigval[ncoef] == 0) stop("Zero eigenvalue of coefficient matrix.")
    logcondition <- log10(eigval[1]) - log10(eigval[ncoef])
    if (logcondition > 12) {
        warning("Near singularity in coefficient matrix.")
        cat(paste("\nEigenvalues range from",eigval[ncoef]," to ",eigval[1]))
    }
}

