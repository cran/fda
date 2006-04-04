plot.basisfd <- function(x, ...) {
  basisobj <- x
	#  plot a basis object
	
	#  Last modified 26 October 2005

    #  check BASISOBJ

    if (!inherits(basisobj, "basisfd")) stop(
        "BASISOBJ is not a basis object.")
		
	dot.args <- list(...)

  	xlabel <- dot.args$xlab
 	if (is.null(xlabel)) xlabel <- ""

  	ylabel <- dot.args$ylab
  	if (is.null(ylabel)) ylabel <- ""

  	cexval <- dot.args$cex
  	if (is.null(cexval)) cexval <- 1

	nbasis   <- basisobj$nbasis
	ltype <- dot.args$lty
	if (is.null(ltype)) ltype=rep(c(1,2,3),nbasis/3)

	nx       <- max(101,10*nbasis)
	rangex   <- basisobj$rangeval
	argvals  <- seq(rangex[1],rangex[2],len=nx)
	basismat <- eval.basis(argvals, basisobj)
	minval   <- min(basismat)
	maxval   <- max(basismat)
	if (minval == maxval) {
    	if (abs(minval) < 1e-1) {
        	minval <- minval - 0.05
        	maxval <- maxval + 0.05
    	} else {
         minval <- minval - 0.05*minval
         maxval <- maxval + 0.05*minval
   		}
	}

	matplot (argvals, basismat, type="l", lty=ltype,
	         xlab=xlabel, ylab=ylabel, cex=cexval,
	         xlim=c(argvals[1],argvals[nx]),
	         ylim=c(minval, maxval))
}

