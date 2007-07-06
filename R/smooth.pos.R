smooth.pos <- function(argvals, y, WfdParobj, wt=rep(1,nobs), conv=1e-4,
                       iterlim=20, dbglev=1) {
# POSFD estimates a positive function fitting a sample of scalar observations.

#  Arguments are:
#  ARGVALS   ... array of function values
#  Y         ... array of argument values
#  WFDPARobj ... functional parameter object defining initial log smooth
#  WT        ... a vector of weights
#  CONV      ... Foldlistergence test value
#  ITERLIM   ... iteration limit for scoring iterations
#  DBGLEV    ... level of output of computation history

#  Returns:
#  WFDOBJ    functional data object defining final smooth function.
#  FLIST      Struct object containing
#               FLIST$f     final log likelihood
#               FLIST$norm  final norm of gradient
#  ITERNUM   Number of iterations
#  ITERHIST  History of iterations

#  last modified 26 October 2005


   if (!(inherits(WfdParobj, "fdPar")))
		stop("Argument WFDPAROBJ not a functional parameter object.")

   lambda <- WfdParobj$lambda
	
	Wfdobj   <- WfdParobj$fd
	Lfdobj   <- WfdParobj$Lfd
	lambda   <- WfdParobj$lambda
	basisobj <- Wfdobj$basis
	nbasis   <- basisobj$nbasis
	rangex   <- basisobj$rangeval

	N  <- length(argvals)
	if (length(y) != N) stop("ARGVALS and Y are not of the same length.")
	
	#  check for argument values out of range
	
	inrng <- (1:N)[argvals >= rangex[1] & argvals <= rangex[2]]
	if (length(inrng) != N)
    	warning("Some values in argvals out of range and not used.")

	argvals <- argvals[inrng]
	y       <- y[inrng]
	nobs    <- length(argvals)

	#  set up some arrays

	climit  <- c(rep(-50,nbasis),rep(400,nbasis))
	cvec0   <- Wfdobj$coefs
	hmat    <- matrix(0,nbasis,nbasis)
	active  <- 1:nbasis
	dbgwrd  <- dbglev > 1

	#  initialize matrix Kmat defining penalty term

	if (lambda > 0)
	  	Kmat <- lambda*eval.penalty(basisobj, Lfdobj)

	#  evaluate log likelihood
	#    and its derivatives with respect to these coefficients

	result <- loglfnpos(argvals, y, basisobj, cvec0, Kmat, wt)
	logl   <- result[[1]]
	Dlogl  <- result[[2]]

	#  compute initial badness of fit measures

	f0    <- -logl
	gvec0 <- -Dlogl
	if (lambda > 0) {
   		gvec0 <- gvec0 +            2*Kmat %*% cvec0
   		f0    <- f0    + t(cvec0) %*% Kmat %*% cvec0
	}
	Foldlist <- list(f = f0, norm = sqrt(mean(gvec0^2)))

	#  compute the initial expected Hessian

	hmat0 <- loglhesspos(argvals, y, basisobj, cvec0, Kmat, wt)
	if (lambda > 0) hmat0 <- hmat0 + 2*Kmat

	#  evaluate the initial update vector for correcting the initial bmat

	deltac   <- -solve(hmat0,gvec0)
	cosangle <- -sum(gvec0*deltac)/sqrt(sum(gvec0^2)*sum(deltac^2))

	#  initialize iteration status arrays

	iternum <- 0
	status <- c(iternum, Foldlist$f, -logl, Foldlist$norm)
	cat("Iteration  Criterion  Neg. Log L  Grad. Norm\n")
	cat("      ")
	cat(format(iternum))
	cat("    ")
	cat(format(status[2:4]))
	cat("\n")
	iterhist <- matrix(0,iterlim+1,length(status))
	iterhist[1,]  <- status
	if (iterlim == 0) {
    	Flist     <- Foldlist
    	iterhist <- iterhist[1,]
		return( list("Wfdobj"=Wfdobj, "Flist"=Flist,
			          "iternum"=iternum, "iterhist"=iterhist) )
	} else {
		gvec <- gvec0
		hmat <- hmat0
	}

	#  -------  Begin iterations  -----------

	STEPMAX <- 5
	MAXSTEP <- 400
	trial   <- 1
	cvec    <- cvec0
	linemat <- matrix(0,3,5)

	for (iter in 1:iterlim) {
   		iternum <- iternum + 1
	   	#  take optimal stepsize
   		dblwrd <- c(FALSE,FALSE)
		limwrd <- FALSE
		stpwrd <- FALSE
		ind    <- 0
	   	#  compute slope
      	Flist <- Foldlist
      	linemat[2,1] <- sum(deltac*gvec)
      	#  normalize search direction vector
      	sdg     <- sqrt(sum(deltac^2))
      	deltac  <- deltac/sdg
      	dgsum   <- sum(deltac)
      	linemat[2,1] <- linemat[2,1]/sdg
      	#  return with stop condition if (initial slope is nonnegative
      	if (linemat[2,1] >= 0) {
        	print("Initial slope nonnegative.")
        	ind <- 3
        	iterhist <- iterhist[1:(iternum+1),]
        	break
      	}
      	#  return successfully if (initial slope is very small
      	if (linemat[2,1] >= -1e-5) {
        	if (dbglev>1) print("Initial slope too small")
        	iterhist <- iterhist[1:(iternum+1),]
        	break
      	}
      	linemat[1,1:4] <- 0
      	linemat[2,1:4] <- linemat[2,1]
      	linemat[3,1:4] <- Foldlist$f
      	stepiter  <- 0
      	if (dbglev > 1) {
			cat("              ")
			cat(format(stepiter))
			cat(format(linemat[,1]))
			cat("\n")
		}
      	ips <- 0
      	#  first step set to trial
      	linemat[1,5]  <- trial
      	#  Main iteration loop for linesrch
      	for (stepiter in 1:STEPMAX) {
        	#  ensure that step does not go beyond limits on parameters
        	limflg  <- 0
        	#  check the step size
        	result <- stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                            climit, active, dbgwrd)
			linemat[1,5] <- result[[1]]
			ind          <- result[[2]]
			limwrd       <- result[[3]]
       	if (linemat[1,5] <= 1e-9) {
          		#  Current step size too small  terminate
          		Flist    <- Foldlist
          		cvecnew <- cvec
          		gvecnew <- gvec
          		if (dbglev > 1) print(paste("Stepsize too small:", linemat[1,5]))
          		if (limflg) ind <- 1 else ind <- 4
          		break
        	}
        	cvecnew <- cvec + linemat[1,5]*deltac
        	#  compute new function value and gradient
			result <- loglfnpos(argvals, y, basisobj, cvecnew, Kmat, wt)
			logl  <- result[[1]]
			Dlogl <- result[[2]]
        	Flist$f  <- -logl
        	gvecnew <- -Dlogl
        	if (lambda > 0) {
            	gvecnew <- gvecnew + 2*Kmat %*% cvecnew
            	Flist$f <- Flist$f + t(cvecnew) %*% Kmat %*% cvecnew
        	}
        	Flist$norm <- sqrt(mean(gvecnew^2))
        	linemat[3,5] <- Flist$f
        	#  compute new directional derivative
        	linemat[2,5] <- sum(deltac*gvecnew)
      		if (dbglev > 1) {
				cat("              ")
				cat(format(stepiter))
				cat(format(linemat[,1]))
				cat("\n")
			}
        	#  compute next step
			result <- stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbgwrd)
			linemat <- result[[1]]
			ips     <- result[[2]]
			ind     <- result[[3]]
			dblwrd  <- result[[4]]
        	trial   <- linemat[1,5]
        	#  ind == 0 implies Foldlistergence
        	if (ind == 0 | ind == 5) break
        	#  end of line search loop
     	}

     	cvec <- cvecnew
     	gvec <- gvecnew
	  	Wfdobj$coefs <- cvec
     	status <- c(iternum, Flist$f, -logl, Flist$norm)
     	iterhist[iter+1,] <- status
		cat("      ")
		cat(format(iternum))
		cat("    ")
		cat(format(status[2:4]))
		cat("\n")
     	#  test for Foldlistergence
     	if (as.numeric(abs(Flist$f - Foldlist$f)) < conv) {
       	iterhist <- iterhist[1:(iternum+1),]
			return( list("Wfdobj"=Wfdobj, "Flist"=Flist,
			             "iternum"=iternum, "iterhist"=iterhist) )
     	}
     	if (Flist$f >= Foldlist$f) break
     	#  compute the Hessian
     	hmat <- loglhesspos(argvals, y, basisobj, cvec, Kmat, wt)
     	if (lambda > 0) hmat <- hmat + 2*Kmat
     	#  evaluate the update vector
     	deltac <- -solve(hmat,gvec)
     	cosangle  <- -sum(gvec*deltac)/sqrt(sum(gvec^2)*sum(deltac^2))
     	if (cosangle < 0) {
       	if (dbglev > 1) print("cos(angle) negative")
       	deltac <- -gvec
     	}
     	Foldlist <- Flist
		#  end of iterations
  	}
	#  compute final normalizing constant
	return( list("Wfdobj"=Wfdobj, "Flist"=Flist,
			      "iternum"=iternum, "iterhist"=iterhist) )
}

#  ---------------------------------------------------------------

loglfnpos <- function(argvals, y, basisobj, cvec, Kmat, wt) {
	#  Computes the log likelihood and its derivative with
	#    respect to the coefficients in CVEC
   	N       <- length(argvals)
   	nbasis  <- basisobj$nbasis
   	phimat  <- getbasismatrix(argvals, basisobj)
	Wvec    <- phimat %*% cvec
	EWvec   <- exp(Wvec)
	res     <- y - EWvec
   	logl    <- -mean(wt*res^2) - t(cvec) %*% Kmat %*% cvec
  	Dlogl   <- 2*crossprod(phimat,wt*res*EWvec)/N - 2*Kmat %*% cvec
	return( list(logl, Dlogl) )
}

#  ---------------------------------------------------------------

loglhesspos <- function(argvals, y, basisobj, cvec, Kmat, wt) {
	#  Computes the expected Hessian
   	N       <- length(argvals)
   	nbasis  <- basisobj$nbasis
   	phimat  <- getbasismatrix(argvals, basisobj)
	Wvec    <- phimat %*% cvec
	EWvec   <- exp(Wvec)
	res     <- y - EWvec
	Dres    <- ((wt*res*EWvec) %*% matrix(1,1,nbasis)) * phimat
  	D2logl  <- 2*crossprod(Dres)/N + 2*Kmat
	return(D2logl)
}
	

