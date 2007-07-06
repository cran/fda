plotfit.fd <- function(y, argvals, fdobj, rng = rangeval,
                       index = 1:nrep, nfine = 101,
                       residual = FALSE, sortwrd = FALSE,
                       titles=NULL,  ylim=NULL, ask=FALSE,
                       type=c("p", "l")[1+residual], xlab=argname,
                       ylab, sub=Sub, col=1:9, lty=1:9, lwd=1,
                       cex.pch=1, ...)
{
#PLOTFIT plots discrete data along with a functional data object for 
#  fitting the data.  It is designed to be used after DATA2FD, 
#  SMOOTH.FD or SMOOTH.BASIS to check the fit of the data offered by
#  the FD object.
  
#  Arguments:
#  Y        ... the data used to generate the fit
#  ARGVALS  ... discrete argument values associated with data
#  fdobj       ... a functional data object for fitting the data
#  RNG      ... a range of argument values to be plotted
#  INDEX    ... an index for plotting subsets of the curves
#       (either sorted or not)
#  NFINE    ... number of points to use for plotting curves
#  RESIDUAL ... if TRUE, the residuals are plotted instead of
#         the data plus curve
#  SORTWRD  ... sort plots by mean square error
#  TITLES   ... vector of title strings for curves

# Last modified 2007.10.03 by Spencer Graves
  
#  Previously modified 20 March 2006

  dots <- list(...)
  if(is.null(titles) && ("main" %in% names(dots)))
    titles <- dots$main
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
#
  casenames <- {
    if(nrep == 1)names(fdnames)[[2]]
    else {
      if(is.null(fdnames[[2]])) dimnames(y)[[2]]
      else fdnames[[2]]
    }
  }
#
  varnames <- {
    if (nvar == 1) names(fdnames)[[3]]
    else {
      if(is.null(fdnames[[3]])) dimnames(y)[[3]]
      else fdnames[[3]]
    }
  }
#  
  if (is.null(argname)) argname <- "Argument Value"
  if (is.null(casenames) || length(casenames) != nrep)
		casenames <- as.character(1:nrep)
  if (is.null( varnames) || length( varnames) != nvar)
		varnames  <- as.character(1:nvar)

#  compute fitted values for evalargs and fine mesh of values
	
  yhat   <- array(eval.fd(argvals, fdobj),c(n,nrep,nvar))
  res    <- y - yhat
  MSE    <- apply(res^2,c(2,3),mean)
  dimnames(MSE) <- list(casenames, varnames)
  MSEsum <- apply(MSE,1,sum)
	
#  compute fitted values for fine mesh of values
	
  xfine <- seq(rng[1], rng[2], len=nfine)
  yfine <- array(eval.fd(xfine, fdobj),c(nfine,nrep,nvar))
	
#  sort cases by MSE if desired
	
  if (sortwrd && nrep > 1) {
    MSEind <- order(MSEsum)
    y      <- y    [,MSEind,, drop=FALSE]
    yhat   <- yhat [,MSEind,, drop=FALSE]
    yfine  <- yfine[,MSEind,, drop=FALSE]
    res    <- res  [,MSEind,, drop=FALSE]
    MSE    <- MSE  [ MSEind,, drop=FALSE]
    casenames  <- casenames[MSEind]
    titles <- titles[MSEind]
#    dim(y)     <- c(n,    nrep,nvar)
#    dim(yhat)  <- c(n,    nrep,nvar)
#    dim(yfine) <- c(nfine,nrep,nvar)
#    dim(res)   <- c(n,    nrep,nvar)
#    dim(MSE)   <- c(      nrep,nvar)
  }
	
#  set up fit and data as 3D arrays, selecting curves in INDEX
	
  y     <- y    [,index,, drop=FALSE]
  yhat  <- yhat [,index,, drop=FALSE]
  res   <- res  [,index,, drop=FALSE]
  yfine <- yfine[,index,, drop=FALSE]
  MSE   <- MSE  [ index,, drop=FALSE]
  casenames <- casenames[index]
  titles <- titles[index]
  nrepi  <- length(index)
#  dim(y)     <- c(n,    nrep,nvar)
#  dim(yhat)  <- c(n,    nrep,nvar)
#  dim(res)   <- c(n,    nrep,nvar)
#  dim(yfine) <- c(nfine,nrep,nvar)
#  dim(MSE)   <- c(      nrep,nvar)
	
# How many plots on a page?   
  nOnOne <- {
    if(ask) min(nrepi, 
      max(length(col), length(lty), length(lwd), length(cex.pch)))
    else nrepi*nvar
  }
# types of plots
  col <- rep(col, length=nOnOne)
  lty <- rep(lty, length=nOnOne)
  lwd <- rep(lwd, length=nOnOne)
  cex.pch <- rep(cex.pch, length=nOnOne)
 
	#  select values in ARGVALS, Y, and YHAT within RNG
	
  argind    <- ((argvals >= rng[1]) & (argvals <= rng[2]))
  argvals   <- argvals[argind]
  casenames <- casenames[argind]
  y    <- y   [argind,,, drop=FALSE]
  yhat <- yhat[argind,,, drop=FALSE]
  res  <- res [argind,,, drop=FALSE]
  n    <- length(argvals)
#	dim(y)    <- c(n,nrep,nvar)
#	dim(yhat) <- c(n,nrep,nvar)
#	dim(res)  <- c(n,nrep,nvar)
	
  xfiind <- ((xfine >= rng[1]) & (xfine <= rng[2]))
  xfine  <- xfine[xfiind]
  yfine  <- yfine[xfiind,,, drop=FALSE]
  nfine  <- length(xfine)
#  dim(yfine) <- c(nfine,nrep,nvar)
		
#  plot the results
	
  ndigit = abs(floor(log10(min(c(MSE)))) - 1)
# 'ask' is controlled by 'nOnOne' ... 
  op <- par(ask=FALSE)
# Don't ask for the first plot,
# but if ask==TRUE, do ask for succeeding plots 
  on.exit(par(op))
# Reset 'ask' on.exit (normal or otherwise)
  iOnOne <- 0
  {
    if (residual) {
#  plot the residuals
#      ylimit <- range(res)
      if(is.null(ylim))ylim <- range(res)
      if(missing(ylab))ylab=rep("Residuals", nrepi)
#      for (i in 1:nrep) {
      for (j in 1:nvar) {
        for (i in 1:nrepi){
#          if (j==1) par(ask = TRUE) else par(ask = FALSE)
          if(iOnOne %in% c(0, nOnOne)[1:(1+ask)]){
#           plot(argvals, res[,i,j], xlim=rng, ylim=ylimit, 
#              xlab=argname, ylab=paste("Residual for",varnames[j]))
            plot(rng, ylim, type="n", xlab=xlab,
                 ylab=ylab[i], ...)
            par(ask=ask)
            abline(h=0, lty=4, lwd=2)
            if(nOnOne==1){
              Sub <- paste("  RMS residual =",
                           round(sqrt(MSE[i,j]),ndigit))
              {
                if (is.null(titles))
                  title(main=casenames[i],sub=sub)
#                 title(main=paste("Case",casenames[i]),
#                   sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
                else title(main=titles[i], sub=sub)
#                 title(main=paste(titles[i]),
#                    sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
              }
            }
            iOnOne <- 0
          }
          iOnOne <- iOnOne+1
          lines(argvals, res[,i,j], type=type,
                col=col[iOnOne], lty=lty[iOnOne],
                lwd=lwd[iOnOne], cex=cex.pch[iOnOne])
        }
      }
    } else {
#  plot the data and fit
#     ylimit <- range(c(c(y),c(yfine)))
      if(is.null(ylim))ylim <- range(c(c(y),c(yfine)))
      if(missing(ylab))ylab <- varnames
#     for (i in 1:nrep) { 
      for (j in 1:nvar) {
        for (i in 1:nrepi) {
#         if (j==1) par(ask = TRUE) else par(ask = FALSE)
          if(iOnOne %in% c(0, nOnOne)[1:(1+ask)]){
#          plot(argvals, y[,i,j], type="p", xlim=rng, ylim=ylimit, col=1,
#               xlab=argname, ylab=varnames[j])
            plot(rng, ylim, type="n", xlab=xlab,
                 ylab=ylab[i], ...)
            par(ask=ask)
            if(nOnOne==1){
              Sub <- paste("  RMS residual =",
                           round(sqrt(MSE[i,j]),ndigit))
              {
                if(is.null(titles)) title(main=casenames[i],
                                          sub=sub)
                else title(main=titles[i], sub=sub)
              }
            }
            iOnOne <- 0
          }
          iOnOne <- iOnOne+1
          points(argvals, y[,i,j], type=type,
                 xlab=xlab, ylab=varnames[j], col=col[iOnOne],
                 lty=lty[iOnOne], lwd=lwd[iOnOne],
                 cex=cex.pch[iOnOne])
          lines(xfine, yfine[,i,j], col=col[iOnOne],
                lty=lty[iOnOne], lwd=lwd[iOnOne])
#          lines(xfine, yfine[,i,j], col=1, lwd=2)
#          if (is.null(titles))
#            title(main=paste("Case",casenames[i]),
#                  sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
#          else 
#            title(main=paste(titles[i]),
#                  sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
        }
      }			
    }
  }
  invisible(NULL)
}
