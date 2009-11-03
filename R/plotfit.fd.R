#plotfit <- function (x, ...){
#  UseMethod("plotfit")
#}

plotfit.fdSmooth <- function(y, argvals, fdSm, rng = NULL,
                       index = NULL, nfine = 101, residual = FALSE,
                       sortwrd = FALSE, titles=NULL,  ylim=NULL,
                       ask=TRUE, type=c("p", "l")[1+residual],
                       xlab=NULL, ylab=NULL, sub=NULL, col=1:9,
                       lty=1:9, lwd=1, cex.pch=1, axes=NULL, ...){
  plotfit.fd(y, argvals, fdSm$fd, rng = rng, index = index,
             nfine = nfine, residual = residual,
             sortwrd = sortwrd, titles=titles,  ylim=ylim,
             ask=ask, type=c("p", "l")[1+residual],
             xlab=xlab, ylab=ylab, sub=sub, col=1:9, lty=1:9,
             lwd=1, cex.pch=1, axes=axes, ...)
}

plotfit.fd <- function(y, argvals, fdobj, rng = NULL,
                       index = NULL, nfine = 101, residual = FALSE,
                       sortwrd = FALSE, titles=NULL,  ylim=NULL,
                       ask=TRUE, type=c("p", "l")[1+residual],
                       xlab=NULL, ylab=NULL, sub=NULL, col=1:9, lty=1:9,
                       lwd=1, cex.pch=1, axes=NULL, ...)
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

# Last modified 2008.12.12 by Spencer Graves
# Previously modified 2008.12.09 by Jim Ramsay
##
## 1.  Basic checks
##
  {
    if(is.null(axes)){
      if(is.null(fdobj$basis$axes)){
        Axes <- TRUE
        axFun <- FALSE
      }
      else{
        if(!inherits(fdobj$basis$axes, 'list'))
          stop('fdobj$basis$axes must be a list;  ',
               'class(fdobj$basis$axes) = ', class(fdobj$basis$axes))
        if(!(inherits(fdobj$basis$axes[[1]], 'character') ||
             inherits(fdobj$basis$axes[[1]], 'function') ) )
          stop('fdobj$basis$axes[[1]] must be either a function or the ',
               'name of a function;  class(fdobj$basis$axes[[1]]) = ',
               class(fdobj$basis$axes[[1]]) )
        Axes <- FALSE
        axFun <- TRUE
        axList <- c(fdobj$basis$axes, ...)
      }
    }
    else{
      if(is.logical(axes)){
        Axes <- axes
        axFun <- FALSE
      }
      else{
        if(!inherits(axes, 'list'))
          stop('axes must be a logical or a list;  class(axes) = ',
               class(axes))
        if(!(inherits(axes[[1]], 'character') ||
             inherits(axes[[1]], 'function') ) )
          stop('axes[[1]] must be either a function or the ',
               'name of a function;  class(axes[[1]]) = ',
               class(axes[[1]]) )
        Axes <- FALSE
        axFun <- TRUE
        axList <- c(axes, ...)
      }
    }
  }
#
  dots <- list(...)
  if(is.null(titles) && ("main" %in% names(dots)))
    titles <- dots$main
  if (!(inherits(fdobj, "fd")))
    stop("Third argument is not a functional data object.")

  basisobj <- fdobj$basis
  if(is.null(rng))rng <- basisobj$rangeval
#
# Are dimensions compatible?
# If yes, extract names ..
  fdnames <- fdobj$fdnames
  yName <- substring(deparse(substitute(y)), 1, 33)
  fdName <- paste(substring(deparse(substitute(fdobj)), 1, 22),
                  "$coef", sep='')
##
## 2.  Use 'checkDims' to reconcile y and fdoj$coef
##
#  The default fdnames may not work well
  defaultNms <- c(fdnames[2], fdnames[3], x='x')
  if((length(defaultNms[[2]])<2) && !is.null(names(defaultNms))
     && !is.na(names(defaultNms)[2]))
    defaultNms[[2]] <- names(defaultNms)[2]
#
  subset <- checkDims3(y, fdobj$coef, defaultNames = defaultNms,
                       xName=yName, yName=fdName)
  y <- subset$x
  fdobj$coef <- subset$y
#
  n <- dim(y)[1]

  argname <- names(fdnames)[[1]]
  if (is.null(argname)) {
      if (is.null(fdnames[[1]])) argname <- "Argument Value"
      else                       argname <- fdnames[[1]]
  }
  if(is.null(xlab))xlab <- argname
  tnames <- dimnames(y)[[1]]
#
  nrep <- dim(y)[2]
  nvar <- dim(y)[3]
#
  curveno <- 1:nrep
#
  casenames <- dimnames(y)[[2]]
  varnames  <- dimnames(y)[[3]]
##
## 2.  Computed fitted values for argvals and a fine mesh
##
  yhat. <- eval.fd(argvals, fdobj)
  yhat   <- as.array3(yhat.)
#  yhat   <- array(eval.fd(argvals, fdobj),c(n,nrep,nvar))
  res    <- y - yhat
  MSE    <- apply(res^2,c(2,3),mean)
  dimnames(MSE) <- list(casenames, varnames)
  MSEsum <- apply(MSE,1,sum)

#  compute fitted values for fine mesh of values

  xfine <- seq(rng[1], rng[2], len=nfine)
  yfine <- array(eval.fd(xfine, fdobj),c(nfine,nrep,nvar))

#  sort cases by MSE if desired?????

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

  if(is.null(index))index <- 1:nrep
#
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
  nOnOne <- 1
#    { if(ask) min(nrepi,
#      max(length(col), length(lty), length(lwd), length(cex.pch)))
#    else nrepi*nvar  }
# types of plots
  if(ask & ((nvar*nrepi/nOnOne) > 1))
    cat('Multiple plots:  Click in the plot to advance to the next')
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
  if(is.null(sub))
    sub <- paste("  RMS residual =", round(sqrt(MSE),ndigit))
  if(length(sub) != length(MSE)){
    warning('length(sub) = ', length(sub), ' != ',
            length(MSE), ' = ', length(MSE), '; forcing equality')
    sub <- rep(sub, length=length(MSE))
  }
  if(is.null(dim(sub))){
#    warning('is.null(dim(sub)); must match dim(MSE) = ',
#            paste(dim(MSE), collapse=', '), ';  forcing equality.')
    dim(sub) <- dim(MSE)
  }
  if(!all(dim(sub)==dim(MSE))){
    warning('dim(sub) = ', dim(sub), " != dim(MSE) = ",
            paste(dim(MSE), collapse=', '), ';  forcing equality.')
    dim(sub) <- dim(MSE)
  }
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
#      if(missing(ylab))ylab=rep("Residuals", nrepi)
      if(is.null(ylab))ylab=paste('Residuals for', varnames)
#      for (i in 1:nrep) {
      for (j in 1:nvar) {
        for (i in 1:nrepi){
#          if (j==1) par(ask = TRUE) else par(ask = FALSE)
          if(iOnOne %in% c(0, nOnOne)[1:(1+ask)]){
#           plot(argvals, res[,i,j], xlim=rng, ylim=ylimit,
#              xlab=argname, ylab=paste("Residual for",varnames[j]))
            if(is.null(xlab))xlab <- argname
            plot(rng, ylim, type="n", xlab=xlab,
                 ylab=ylab[j], axes=Axes, ...)
            if(axFun)
              do.call(axList[[1]], axList[-1])
#            axis(1)
#            axis(2)
            par(ask=ask)
            abline(h=0, lty=4, lwd=2)
            if(nOnOne==1){
              {
                if (is.null(titles))
                  title(main=casenames[i],sub=sub[i, j])
#                 title(main=paste("Case",casenames[i]),
#                   sub =paste("  RMS residual =",round(sqrt(MSE[i,j]),ndigit)))
                else title(main=titles[i], sub=sub[i, j])
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
      varname <- names(fdnames)[[3]]
      if (is.null(varname)) {
          if (is.null(fdnames[[3]])) rep("Function Value", nvar)
          else                       varname <- fdnames[[3]]
      }
      if(is.null(ylab))ylab <- varname
#     for (i in 1:nrep) {
      for (j in 1:nvar) {
        for (i in 1:nrepi) {
#         if (j==1) par(ask = TRUE) else par(ask = FALSE)
          if(iOnOne %in% c(0, nOnOne)[1:(1+ask)]){
#          plot(argvals, y[,i,j], type="p", xlim=rng, ylim=ylimit, col=1,
#               xlab=argname, ylab=varnames[j])
            plot(rng, ylim, type="n", xlab=xlab,
                 ylab=ylab[j], axes=Axes, ...)
            if(axFun)
              do.call(axList[[1]], axList[-1])
#            axis(1)
#            axis(2)
            par(ask=ask)
            if(nOnOne==1){
              {
                if(is.null(titles)) title(main=casenames[i],
                                          sub=sub[i, j])
                else title(main=titles[i], sub=sub[i, j])
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
