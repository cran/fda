plotfit.fd <- function(x, argvals, fdobj, rng = rangeval, index = 1:nrep, nfine = 101,
                       residual = FALSE, sortwrd = FALSE, titles=NULL, ...)
{
  y <- x
	#PLOTFIT plots discrete data along with a functional data object for fitting the
	#  data.  It is designed to be used after DATA2FD, SMOOTH.FD or SMOOTH.BASIS to
	#  check the fit of the data offered by the FD object.
	#  Arguments:
	#  Y        ... the data used to generate the fit
	#  ARGVALS  ... discrete argument values associated with data
	#  fdobj       ... a functional data object for fitting the data
	#  RNG      ... a range of argument values to be plotted
	#  INDEX    ... an index for plotting subsets of the curves (either sorted or not)
	#  NFINE    ... number of points to use for plotting curves
	#  RESIDUAL ... if TRUE, the residuals are plotted instead of the data plus curve
	#  SORTWRD  ... sort plots by mean square error
	#  TITLES   ... vector of title strings for curves
	
	#  Last modified 20 March 2006
	
   if (!(inherits(fdobj, "fd"))) stop(
		"Third argument is not a functional data object.")

	basisobj <- fdobj$basis
	rangeval <- basisobj$rangeval
	
	coef  <- fdobj$coefs
	coefd <- dim(coef)
	ndim  <- length(coefd)
	
	y <- as.array(y)
	n <- dim(y)[1]
		
	dimnames(y) <- NULL
	if (ndim < 2) nrep <- 1 else nrep <- coefd[2]
	if (ndim < 3) nvar <- 1 else nvar <- coefd[3]
	dim(y) <- c(n, nrep, nvar)
	
	curveno <- 1:nrep
	
	fdnames <- fdobj$fdnames
	argname <- names(fdnames)[[1]]
	if (nrep == 1) casenames <- names(fdnames)[[2]] else casenames <- fdnames[[2]]
	if (nvar == 1) varnames  <- names(fdnames)[[3]] else varnames  <- fdnames[[3]]
	if (is.null(argname)) argname <- "Argument Value"
	if (is.null(casenames) || length(casenames) != nrep)
		casenames <- as.character(1:nrep)
	if (is.null( varnames) || length( varnames) != nvar)
		varnames  <- as.character(1:nvar)

	#  compute fitted values for evalargs and fine mesh of values
	
	yhat   <- array(eval.fd(argvals, fdobj),c(n,nrep,nvar))
	res    <- y - yhat
	MSE    <- apply(res^2,c(2,3),mean)
	MSEsum <- apply(MSE,1,sum)
	
	#  compute fitted values for fine mesh of values
	
	xfine <- seq(rng[1], rng[2], len=nfine)
	yfine <- array(eval.fd(xfine, fdobj),c(nfine,nrep,nvar))
	
	#  sort cases by MSE if desired
	
	if (sortwrd && nrep > 1) {
		MSEind <- order(MSEsum)
		y      <- y    [,MSEind,]
		yhat   <- yhat [,MSEind,]
		yfine  <- yfine[,MSEind,]
		res    <- res  [,MSEind,]
		MSE    <- MSE  [ MSEind,]
		casenames  <- casenames[MSEind]
		dim(y)     <- c(n,    nrep,nvar)
		dim(yhat)  <- c(n,    nrep,nvar)
		dim(yfine) <- c(nfine,nrep,nvar)
		dim(res)   <- c(n,    nrep,nvar)
		dim(MSE)   <- c(      nrep,nvar)
	}
	
	#  set up fit and data as 3D arrays, selecting curves in INDEX
	
	y     <- y    [,index,]
	yhat  <- yhat [,index,]
	res   <- res  [,index,]
	yfine <- yfine[,index,]
	MSE   <- MSE  [ index,]
	casenames <- casenames[index]
	nrep  <- length(index)
	dim(y)     <- c(n,    nrep,nvar)
	dim(yhat)  <- c(n,    nrep,nvar)
	dim(res)   <- c(n,    nrep,nvar)
	dim(yfine) <- c(nfine,nrep,nvar)
	dim(MSE)   <- c(      nrep,nvar)
	
	#  select values in ARGVALS, Y, and YHAT within RNG
	
	argind    <- argvals >= rng[1] & argvals <= rng[2]
	argvals   <- argvals[argind]
	casenames <- casenames[argind]
	y    <- y   [argind,,]
	yhat <- yhat[argind,,]
	res  <- res [argind,,]
	n    <- length(argvals)
	dim(y)    <- c(n,nrep,nvar)
	dim(yhat) <- c(n,nrep,nvar)
	dim(res)  <- c(n,nrep,nvar)
	
	xfiind <- xfine >= rng[1] & xfine <= rng[2]
	xfine  <- xfine[xfiind]
	yfine  <- yfine[xfiind,,]
	nfine  <- length(xfine)
	dim(yfine) <- c(nfine,nrep,nvar)
		
	#  plot the results
	
	ndigit = abs(floor(log10(min(c(MSE)))) - 1)
	if (residual) {
	    #  plot the residuals
	    ylimit <- range(res)
	    for (i in 1:nrep) {
              for (j in 1:nvar) {
                #if (j==1) par(ask = TRUE) else par(ask = FALSE)
		    plot(argvals, res[,i,j], xlim=rng, ylim=ylimit, 
		          xlab=argname, ylab=paste("Residual for",varnames[j]))
		    abline(h=0, lty=4, lwd=2)
		    if (is.null(titles))
		    	  title(main=paste("Case",casenames[i]),
		              sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
		    else title(main=paste(titles[i]),
		               sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
              }
	    }			
	} else {
	    #  plot the data and fit
	    ylimit <- range(c(c(y),c(yfine)))
	    for (i in 1:nrep) { 
              for (j in 1:nvar) {
                #if (j==1) par(ask = TRUE) else par(ask = FALSE)
		    plot(argvals, y[,i,j], type="p", xlim=rng, ylim=ylimit, col=1,
		          xlab=argname, ylab=varnames[j])
		    lines(xfine, yfine[,i,j], col=1, lwd=2)
		    if (is.null(titles))
		    	title(main=paste("Case",casenames[i]),
		            sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
		    else 
                  title(main=paste(titles[i]),
		            sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
              }
	    }			
	}
	invisible(NULL)
}
