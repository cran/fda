plot.fdSmooth <- function(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                          xlim=NULL, ylim=NULL, xlab=NULL,
                          ylab=NULL, ask=FALSE, nx=201, ...){
  plot(x$fd, y, Lfdobj=Lfdobj, href=href, titles=titles,
              xlim=xlim, ylim=ylim, xlab=xlab,
              ylab=ylab, ask=ask, nx=nx, ...)
}

plot.fd <- function(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL,
                    ylab=NULL, ask=FALSE, nx=201, ...)
{
#  -----------------------------------------------------------------------
#       plot for fd class
#  -----------------------------------------------------------------------
#  Plot a functional data object fdobj.
#  Arguments:
#  x = fdobj     ... a functional data object
#  Lfdobj    ... linear differental operator to be applied to fdobj before
#             plotting
#  HREF   ... If TRUE, a horizontal dotted line through 0 is plotted.
# ...   
#  ASK, if TRUE, causes the curves to be displayed one at a time.
#  NX     ... The number of sampling points to use for
  #             plotting.  (default 101)
  
  #  The remaining optional arguments are the same as those available
  #     in the regular "plot" function.

  #  Note that for multivariate fdobj, a suitable matrix of plots
  #    must be set up before calling plot by using something such as
  #    par(mfrow=c(1,nvar),pty="s")

  # last modified 2008.06.28 by Spencer Graves 
  #  previously modified 2007.05.03 and 6 March 2007

  # fdobj, Lfdobj=0, href=TRUE, nx=201, titles=NULL,
  #                   xlab=xlabel, ylab=ylabel,
  #                   xlim=rangex, ylim=rangey, ask=FALSE

  #  check fdobj

  fdobj <- x
  if (!(inherits(fdobj, "fd"))) stop(
		"First argument is not a functional data object.")

  #  check Lfdobj

#  if (missing(Lfdobj)) Lfdobj <- 0
  Lfdobj <- int2Lfd(Lfdobj)
  if (!inherits(Lfdobj, "Lfd")) stop(
      "Second argument is not a linear differential operator.")

  #  check href
# if (missing(href)) href <- TRUE

  #  check ask
# if (missing(ask)) ask <- FALSE

  #  extract dimension information

  coef   <- fdobj$coefs
  coefd  <- dim(coef)
  ndim   <- length(coefd)
# Number of basis functions   
  nbasis <- coefd[1]
# Number of functional observations   
  nrep   <- coefd[2]
  if (ndim > 2) nvar <- coefd[3] else nvar <- 1

  #  get basis information

  basisobj <- fdobj$basis
  rangex   <- basisobj$rangeval
  if(is.null(xlim))xlim <- rangex 
  #  check xlim
# if (missing(xlim)) xlim <- rangex

  #  set up a set of argument values for the plot

  if (missing(y)) {
#    y <- 201
    y <- nx
  } else {
    y <- as.vector(y)
  }

  if (length(y) == 1) {
    if (y >= 1) { 
      y <- seq(rangex[1],rangex[2],len=floor(y))
    } else {
      stop("'y' a single number less than one.")
    }
  }
  if (min(y) < rangex[1] || max(y) > rangex[2]) stop(
    "Values in Y are outside the basis range.")

  #  evaluate LFDOBJ(FDOBJ) at the argument values

  fdmat    <- eval.fd(y, fdobj, Lfdobj)
  rangey   <- range(c(fdmat))
  if(is.null(ylim))ylim <- rangey
  #  check ylim
# if (missing(ylim)) ylim <- rangey

  xlabel   <- names(fdobj$fdnames)[[1]]
  ylabel   <- names(fdobj$fdnames)[[3]]
  if (is.character(xlabel) == FALSE) xlabel <- ""
  if (is.character(ylabel) == FALSE) ylabel <- ""

  #  check xlab and ylab
  if(is.null(xlab))xlab <- xlabel
  if(is.null(ylab))ylab <- ylabel
#  if (missing(xlab)) xlab <- xlabel
#  if (missing(ylab)) ylab <- ylabel
  crvnames <- fdobj$fdnames[[2]]
  varnames <- fdobj$fdnames[[3]]
# Don't ask for the first plot, but do for later plots if(ask)   
  op <- par(ask=FALSE)
# Don't ask for the first plot,
# but if ask==TRUE, do ask for succeeding plots 
  on.exit(par(op))
# A single line?  
  if (ndim < 2) {
    plot (y, fdmat, type="l", xlim=xlim, ylim=ylim,
          xlab=xlab, ylab=ylab, ...)
#    if (zerofind(fdmat) && href) abline(h=0,lty=2)
    if (zerofind(ylim) && href) abline(h=0,lty=2)
  }
# Several copies of one function?    
  if (ndim ==2 ) {
    if (!ask) {
      matplot(y, fdmat, type="l", xlim=xlim, ylim=ylim,
           		xlab=xlab, ylab=ylab, ...)
#      if (zerofind(fdmat) && href) abline(h=0,lty=2)
      if (zerofind(ylim) && href) abline(h=0,lty=2)
    } else  {
      cat('Multiple plots:  Click in the plot to advance to the next') 
#      op <- par(ask = TRUE)
#      on.exit(par(op))
      for (irep in 1:nrep) {
        plot (y, fdmat[,irep], type="l", xlim=xlim, ylim=ylim,
                xlab=xlab, ylab=ylab, ...)
        par(ask=ask)
#        if (zerofind(fdmat[,irep]) && href) abline(h=0,lty=2)
        if (zerofind(ylim) && href) abline(h=0,lty=2)
        if (!is.null(titles)) title(titles[irep])
        else title(paste(crvnames[irep]))
#        else title(paste("Curve",irep,crvnames[irep]), line=0.2)
#       ... "line=0.2" to allow "main" in "..."         
#        mtext("Click in graph to see next plot", side=3, outer=FALSE)
#        text("",locator(1))
      }
    }
  }
# Possibly multiple copies of different functions   
  if (ndim == 3) {
    if (!ask) {
      for (ivar in 1:nvar) {
        matplot (y, fdmat[,,ivar], type="l", xlim=xlim, ylim=ylim,
                 xlab=xlab, ylab=ylab,
                 main=varnames[ivar], ask=FALSE, ...)
#        if (zerofind(fdmat[,,ivar]) && href) abline(h=0,lty=2)
        if (zerofind(ylim) && href) abline(h=0,lty=2)
      }
    } else {
      for (irep in 1:nrep) {
#            if (ivar==1){op <- par(ask = TRUE); on.exit(par(op)) }
#            else { op <- par(ask = FALSE); on.exit(par(op)) }
        for (ivar in 1:nvar) {
          plot(y,fdmat[,irep,ivar],type="l", xlim=xlim, ylim=ylim,
                xlab=xlab, ylab=ylab, ...)
          par(ask=ask)
#          if (zerofind(fdmat[,irep,ivar]) && href) abline(h=0,lty=2)
          if (zerofind(ylim) && href) abline(h=0,lty=2)
          if (!is.null(titles)) title(titles[irep])
          else title(paste("Curve", irep, varnames[ivar]))
#          mtext("Click in graph to see next plot", side=3, outer=FALSE)
#          text("",locator(1))
        }
      }
    }
  }
#  invisible(NULL)
# This used to return 'invisible(NULL)'.
# However, with R 2.7.0 under XEmacs with ESS,
# it sometimes failed to plot.  I changed the return value,
# and the problem disappeared.  
  'done'
}

zerofind <- function(fmat)
{
  frng <- range(fmat)
  if (frng[1] <= 0 && frng[2] >= 0) zeroin <- TRUE else zeroin <- FALSE
  return(zeroin)
}
