plotFd <- function(fd, Lfd=0, matplt=TRUE, href=TRUE, nx=101,
                    xlab=xlabel, ylab=ylabel, 
                    xlim=rangex, ylim=rangey, ...)
{
  #  Plot a functional data object FD.
  #  Arguments:
  #  FD     ... a functional data object
  #  LFD    ... linear differental operator to be applied to FD before
  #             plotting
  #  MATPLT ... If T, all curves are plotted in a single plot.
  #             Otherwise, each curve is plotted separately, and the
  #             next curve is plotted when the mouse is clicked.
  #  HREF   ... If T, a horizontal dotted line through 0 is plotted.
  #  NX     ... The number of sampling points to use for
  #             plotting.  (default 101)
  #  The remaining optional arguments are the same as those available
  #     in the regular 'plot' function.

  #  Note that for multivariate FD, a suitable matrix of plots
  #    must be set up before calling plot by using something such as
  #    par(mfrow=c(1,nvar),pty='s')

  #  Last modified 4 July 2001

  if (!(inherits(fd, "fd"))) stop('First argument is not a functional data object.')
  if (!is.Lfd(Lfd)) stop(
      "Second argument is not a linear differential operator.")

  coef   <- getcoef(fd)
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  nbasis <- coefd[1]
  nrep   <- coefd[2]
  if (ndim > 2) nvar <- coefd[3] else nvar <- 1

  basisfd <- getbasis(fd)
  rangex  <- basisfd$rangeval
  x       <- seq(rangex[1],rangex[2],length=nx)
  fdmat   <- eval.fd(x, fd, Lfd)
  rangey  <- range(c(fdmat))

  xlabel   <- names(fd$fdnames)[[1]]
  ylabel   <- names(fd$fdnames)[[3]]
  if (is.character(xlabel) == FALSE) xlabel <- ''
  if (is.character(ylabel) == FALSE) ylabel <- ''
  crvnames <- fd$fdnames[[2]]
  varnames <- fd$fdnames[[3]]

  if (ndim < 2) {
    plot (x, fdmat, type='l', xlim=xlim, ylim=ylim, 
          xlab=xlab, ylab=ylab, ...)
    if (zerofind(fdmat) && href) abline(h=0,lty=2)
  }
  if (ndim ==2 ) {
    if (matplt) {
      matplot (x, fdmat, type='l', xlim=xlim, ylim=ylim, 
           xlab=xlab, ylab=ylab, ...)
      if (zerofind(fdmat) && href) abline(h=0,lty=2)
    } else  {
      for (irep in 1:nrep) {
        plot (x, fdmat[,irep], type='l', xlim=xlim, ylim=ylim, 
                xlab=xlab, ylab=ylab,
                main=paste('Curve',irep,crvnames[irep]),...)
        if (zerofind(fdmat[,irep]) && href) abline(h=0,lty=2)
        mtext('Click to advance to next set of plots',side=3,line=-3,outer=TRUE)
        text(locator(1),'')
      }
    }
  }
  if (ndim == 3) {
    if (matplt) {
      for (ivar in 1:nvar) {
        matplot (x, fdmat[,,ivar], type='l', xlim=xlim, ylim=ylim, 
                 xlab=xlab, ylab=ylab,
                 main=varnames[ivar],...)
        if (zerofind(fdmat[,,ivar]) && href) abline(h=0,lty=2)
      }
    }
    if (!matplt)  {
      for (irep in 1:nrep) {
        for (ivar in 1:nvar) {
          plot (x,fdmat[,irep,ivar],type='l', xlim=xlim, ylim=ylim, 
                xlab=xlab, ylab=ylab,
                main=paste('Curve', irep, varnames[ivar]),...)
          if (zerofind(fdmat[,irep,ivar]) && href) abline(h=0,lty=2)
        }
        mtext('Click to advance to next set of plots',side=3,line=-3,outer=TRUE)
        text(locator(1),'')
      }
    }
  }
  invisible()
}

zerofind <- function(fmat)
{
  frng <- range(fmat)
  if (frng[1] <= 0 && frng[2] >= 0) zeroin <- TRUE else zeroin <- FALSE
  return(zeroin)
}
