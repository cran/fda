fRegress <- function(yfdPar, xfdlist, betalist) {

#function [betaestlist, yhatfdobj, betastderrlist, bvar, c2bmap] =
#    fRegress(yfdPar, xfdlist, betalist, y2cmap, sigmae)
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
#  BETALIST ... a list object of length p with each list
#               containing a functional parameter object for
#               the corresponding regression function.
#
#  Returns:
#  FREGRESSLIST   ...  A list containing six members with names:
#    yfdPar      ... first  argument of FREGRESS
#    xfdlist     ... second argument of FREGRESS
#    betalist    ... third  argument of FREGRESS
#    betaestlist ... estimated regression functions
#    yhatfdobj   ... functional data object containing fitted functions
#    Cmatinv     ... inverse of the coefficient matrix, needed for
#                    function FREGRESS.STDERR that computes standard errors

#  Last modified 30 October 2005

#  get number of independent variables

p <- length(xfdlist)

#  check BETALIST

if (length(betalist) != p)  {
	cat(paste("\nNumber of regression coefficients does not match\n",
		       "number of independent variables."))
	stop("")
}

for (j in 1:p) {
	betafdParj <- betalist[[j]]
	if (inherits(betafdParj, "fd")) {
		betafdParj    <- fdPar(betafdParj)
		betalist[[j]] <- betafdParj
	}
	if (!inherits(betafdParj, "fdPar")) stop(
		paste("BETALIST[[",j,"]] is not a FDPAR object."))
}


if (inherits(yfdPar, "fdPar") || inherits(yfdPar, "fd")) {
	
    #  ----------------------------------------------------------------
    #           YFDPAR is functional for a functional parameter
    #  ----------------------------------------------------------------

    if (inherits(yfdPar, "fd")) yfdPar <- fdPar(yfdPar)

    yfdobj  <- yfdPar$fd
    ylambda <- yfdPar$lambda
    yLfdobj <- yfdPar$Lfd
    ycoef   <- yfdobj$coefs
    if (length(dim(ycoef)) > 2) stop("YFDOBJ from YFDPAR is not univariate.")
    N         <- dim(ycoef)[2]
    ybasisobj <- yfdobj$basis
    rangeval  <- ybasisobj$rangeval
    ynbasis   <- ybasisobj$nbasis
    nfine     <- max(501,10*ynbasis+1)
    tfine     <- seq(rangeval[1], rangeval[2], len=nfine)
    ywtvec    <- matrix(1,nfine,nfine)
    deltat    <- tfine[2] - tfine[1]
    ymat      <- eval.fd(tfine, yfdobj)

    #  check each xfdlist member.  If the object is a vector of length N,
    #  it is converted to a functional data object with a
    #  constant basis

    onebasis <- create.constant.basis(rangeval)

    for (j in 1:p) {
        xfdj <- xfdlist[[j]]
        if (inherits(xfdj,"fd")) {
            xcoef <- xfdj$coefs
            if (length(dim(xcoef)) > 2) stop(
					"Covariate is not univariate.")
            rangevalx  <- xfdj$basis$rangeval
            if (any(rangevalx != rangeval))  stop(
					"Range for covariate does not match that of YFDOBJ")
        }
        else if (inherits(xfdj, "numeric"))
				#  convert to a functional data object with cosntant basis
					xfdlist[[j]] <- fd(matrix(xfdj,1,N), onebasis)
        else  stop(paste("Covariate", j, "is neither a functional nor",
                         " a multivariate object."))

        #  check size of coefficient array
        xfdj  <- xfdlist[[j]]
        coefj <- xfdj$coefs
        Nj <- dim(coefj)[2]
        if (Nj != N) stop("Incorrect number of replications in XFDLIST")
    }

    #  set up a matrix of values of covariates over a fine mesh

    xmat <- array(0,c(nfine, N, p))

    betamatlist <- vector("list",p)
    for (j in 1:p) {
        xmatj            <- eval.fd(tfine, xfdlist[[j]])
        xmat[,,j]        <- xmatj
        betafdParj       <- betalist[[j]]
        betabasisj       <- betafdParj$fd$basis
        betamatlist[[j]] <- eval.basis(tfine, betabasisj)
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
            ncoefj     <- betafdParj$fd$basis$nbasis
            #  row indices of CMAT and DMAT to fill
            mj1    <- mj2 + 1
            mj2    <- mj2 + ncoefj
            indexj <- mj1:mj2
            #  compute right side of equation DMAT
        	  xywtvec  <- apply(xmat[,,j]*ymat,1,sum)
            betamatj <- betamatlist[[j]]
            Dmatj    <- deltat*apply(betamatj*xywtvec,2,sum)
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
                    if (N > 1)
            		       xxwtvec  <- apply(xmat[,,j]*xmat[,,k],1,sum)
                    else
                        xxwtvec <- xmat[,1,j]*xmat[,1,k]
                    betamatk <- betamatlist[[k]]
                    Cmatjk   <- deltat*crossprod(betamatj*xxwtvec,betamatk)
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

    #  check Cmat for singularity

    eigchk(Cmat)

    #  solve for coefficients defining BETA

    Cmat    <- (Cmat+t(Cmat))/2
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
        betafdj$coefs <- coefj
        betafdParj$fd <- betafdj
        }
        betaestlist[[j]] <- betafdParj
    }

    #  set up fd objects for predicted values in YHATFDOBJ

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

 }
 else if (inherits(yfdPar,"numeric")) {

    #  ----------------------------------------------------------------
    #                   YFDPAR is scalar or multivariate
    #  ----------------------------------------------------------------

    ymat <- as.matrix(yfdPar)
    N    <- dim(ymat)[1]

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
            Zmat  <- cbind(Zmat,Zmatj)
            ncoefj <- dim(Zmatj)[2]
            Rmatj <- matrix(0,ncoefj,ncoefj)
            Rmat  <- rbind(cbind(Rmat, matrix(0,ncoef,ncoefj)),
                           cbind(matrix(0,ncoefj,ncoef), Rmatj))
            ncoef <- ncoef + ncoefj
          }
        else stop(
			paste("XFDLIST[[",j,
			      "is neither a functional nor a multivariate object."))

    }

    #  -----------------------------------------------------------
    #          set up the linear equations for the solution
    #  -----------------------------------------------------------

    #  solve for coefficients defining BETA

    Cmat <- t(Zmat) %*% Zmat + Rmat
    Dmat <- t(Zmat) %*% ymat

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
            betaestfdj$coefs <- betacoefj
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

 }
 else
    #  YFDOBJ is neither functional nor multivariate
    stop("YFDOBJ is neither functional nor multivariate.")

fRegressList <-
       list(yfdPar      = yfdPar,
            xfdlist     = xfdlist,
            betalist    = betalist,
            betaestlist = betaestlist,
            yhatfdobj   = yhatfdobj,
            Cmatinv     = Cmatinv)

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

