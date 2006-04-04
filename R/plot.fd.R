plot.fd <- function(x, Lfdobj=0, href=TRUE, nx=201, titles=NULL,
                    xlab=xlabel, ylab=ylabel,
                    xlim=rangex, ylim=rangey, ask=FALSE, ...)
{
  fdobj <- x
  #  Plot a functional data object fdobj.
  #  Arguments:
  #  fdobj     ... a functional data object
  #  Lfdobj    ... linear differental operator to be applied to fdobj before
  #             plotting
  #  HREF   ... If TRUE, a horizontal dotted line through 0 is plotted.
  #  NX     ... The number of sampling points to use for
  #             plotting.  (default 101)
  #  ASK    ... If True, each replication is plotted in turn, with 
  #             the request "Waiting to confirm page change"
  #             appearing in the R Console window.  A mouse click
  #             brings up the next plot.  Default is False, which
  #             displays all the curves at once.

  #  The remaining optional arguments are the same as those available
  #     in the regular "plot" function.

  #  The argument ASK, if TRUE, causes the curves to be displayed one at a time.

  #  Note that for multivariate fdobj, a suitable matrix of plots
  #    must be set up before calling plot by using something such as
  #    par(mfrow=c(1,nvar),pty="s")

  #  Last modified 21 March 2006

  #  check fdobj

  if (!(inherits(fdobj, "fd"))) stop(
		"First argument is not a functional data object.")

  #  check Lfdobj

  Lfdobj = int2Lfd(Lfdobj)

  if (!inherits(Lfdobj, "Lfd")) stop(
      "Second argument is not a linear differential operator.")

  coef   <- fdobj$coefs
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  nbasis <- coefd[1]
  nrep   <- coefd[2]
  if (ndim > 2) nvar <- coefd[3] else nvar <- 1

  basisobj <- fdobj$basis
  rangex   <- basisobj$rangeval
  x        <- seq(rangex[1],rangex[2],length=nx)
  fdmat    <- eval.fd(x, fdobj, Lfdobj)
  rangey   <- range(c(fdmat))

  xlabel   <- names(fdobj$fdnames)[[1]]
  ylabel   <- names(fdobj$fdnames)[[3]]
  if (is.character(xlabel) == FALSE) xlabel <- ""
  if (is.character(ylabel) == FALSE) ylabel <- ""
  crvnames <- fdobj$fdnames[[2]]
  varnames <- fdobj$fdnames[[3]]

  if (ndim < 2) {
    plot (x, fdmat, type="l", xlim=xlim, ylim=ylim,
          xlab=xlab, ylab=ylab, ...)
    if (zerofind(fdmat) && href) abline(h=0,lty=2)
  }
  if (ndim ==2 ) {
    if (!ask) {
      matplot(x, fdmat, type="l", xlim=xlim, ylim=ylim,
           		xlab=xlab, ylab=ylab, ...)
      if (zerofind(fdmat) && href) abline(h=0,lty=2)
    } else  {
      #par(ask = TRUE)
      for (irep in 1:nrep) {
        plot (x, fdmat[,irep], type="l", xlim=xlim, ylim=ylim,
                xlab=xlab, ylab=ylab, ...)
        if (zerofind(fdmat[,irep]) && href) abline(h=0,lty=2)
        if (!is.null(titles)) title(titles[irep])
        else title(paste("Curve",irep,crvnames[irep]))
      }
    }
  }
  if (ndim == 3) {
    if (!ask) {
      for (ivar in 1:nvar) {
        matplot (x, fdmat[,,ivar], type="l", xlim=xlim, ylim=ylim,
                 xlab=xlab, ylab=ylab,
                 main=varnames[ivar], ...)
        if (zerofind(fdmat[,,ivar]) && href) abline(h=0,lty=2)
      }
    } else {
      for (irep in 1:nrep) {
        for (ivar in 1:nvar) {
          #if (ivar==1) par(ask = TRUE) else par(ask = FALSE)
          plot(x,fdmat[,irep,ivar],type="l", xlim=xlim, ylim=ylim,
                xlab=xlab, ylab=ylab, ...)
          if (zerofind(fdmat[,irep,ivar]) && href) abline(h=0,lty=2)
          if (!is.null(titles)) title(titles[irep])
          else title(paste("Curve", irep, varnames[ivar]))
        }
      }
    }
  }
  invisible(NULL)
}

zerofind <- function(fmat)
{
  frng <- range(fmat)
  if (frng[1] <= 0 && frng[2] >= 0) zeroin <- TRUE else zeroin <- FALSE
  return(zeroin)
}

