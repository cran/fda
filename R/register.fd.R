register.fd <- function(y0fd, yfd, WfdParobj,
                        conv=1e-4, iterlim=20, dbglev=1, periodic=FALSE, crit=2)
{
#REGISTERFD registers a set of curves YFD to a target function Y0FD.
#  Arguments are:
#  Y0FD      ... Functional data object for target function.  It must
#                contain a single curve, but this single curve
#                can be multivariate.
#  YFD       ... Functional data object for functions to be registered
#  WFDPAROBJ ... Functional parameter object for function W defining warping fns
#                Its coefficients are the starting values used in the
#                iterative computation of the final warping fns.
#                NB:  The first coefficient is is NOT used.
#                For both B-spline and Fourier bases, this first
#                coefficient determines the constant term in the expansion,
#                and, since a register function is normalized, this term
#                is, in effect, eliminated or has no influence on the
#                result.  This first position is used, however, to
#                contain the shift parameter in case the data are
#                treated as periodic.  At the end of the calculations,
#                the shift parameter is returned separately.
#  CONV    ... Convergence criterion
#  ITERLIM  .. iteration limit for scoring iterations
#  DBGLEV  ... Level of output of computation history
#  PERIODIC .. If one, curves are periodic and a shift parameter is fit.
#              Initial value for shift parameter is taken to be 0.
#              The periodic option should ONLY be used with a Fourier
#              basis for the target function Y0FD, the functions to be
#              registered, YFD, and the functions WFD0 defining the
#              time-warping functions.
#  CRIT    ... If 1 least squares, if 2 log eigenvalue ratio.  Default is 1.
#                Default:  0
#  Returns:
#  REGSTR  ...  A list with fields
#    REGSTR$REGFD ... A functional data object for the registered curves
#    REGSTR$WFD   ... A Functional data object for function W defining
#                         warping fns
#    REGSTR$SHIFT ... Shift parameter value if curves are periodic

#  last modified 13 March 2007

#  check classes of first two arguments

if (!(inherits(y0fd, "fd")))
	stop("First argument is not a functional data object.")

if (!(inherits(yfd, "fd")))
	stop("Second argument is not a functional data object.")
	
#  Check target function to determine number of variables.

#  If the coefficient matrix is 3-dimensional, then the second
#  dimension size must be 1, and if not, an error will be declared.
#  If it is 1, however, then it the size of the third dimension will
#  be taken to be the number of variables for a multivariate set
#  of curves.  A 3-dimensional coefficient matrix for the target 
#  will then be contracted to two dimensions.  

y0coefs <- y0fd$coefs
y0dim   <- dim(y0coefs)
ndimy0  <- length(y0dim)
if (ndimy0 == 3) {
   if (y0dim[2] > 1) {
       stop("Y0FD is not a single function.")
   } else {
       y0coefs <- y0coefs[,1,]
       y0fd$coefs <- y0coefs
   }
}

#  The target function coefficient matrix is now taken to be 
#  2-dimensional, and the size of the second dimension is taken to be
#  the number of variables to be registered. If this size
#  is greater than one, the functional data objects to be registered 
#  are assumed to be multivariate.

y0coefs <- y0fd$coefs
y0dim   <- dim(y0coefs)
nvar    <- y0dim[2]
    
#  check functions to be registered

ydim   <- dim(yfd$coefs)
ncurve <- ydim[2]
ndimy  <- length(ydim)

if (ndimy > 3) stop("YFD is more than 3-dimensional.")

if (ndimy == 2) {
   if (nvar > 1) stop(
       paste("Y0FD indicates function are multivariate,",
             "but is YFD is only 2-dimensional."))
}

#  Extract basis information from YFD

ybasis  <- yfd$basis
ynbasis <- ybasis$nbasis
if (periodic && !(ybasis$type == "fourier"))
      stop("PERIODIC is true, basis not fourier type.")

#  check functions W defining warping functions

Wfd0   <- WfdParobj$fd
wcoef  <- Wfd0$coefs
wbasis <- Wfd0$basis
nbasis <- wbasis$nbasis
wtype  <- wbasis$type
rangex <- wbasis$rangeval
wdim   <- dim(wcoef)
ncoef  <- wdim[1]
ndimw  <- length(wdim)
if (wdim[ndimw] == 1) ndimw <- ndimw - 1
if (ndimw == 1 && ncurve > 1)
      stop("WFD and YFD do not have the same dimensions.")
if (ndimw == 2 && wdim[2] != ncurve)
      stop("WFD and YFD do not have the same dimensions.")
if (ndimw > 2)  stop("WFD is not univariate.")

#  set up a fine mesh of argument values

NFINEMIN <- 101
nfine <- 10*ynbasis + 1
if (nfine < NFINEMIN) nfine <- NFINEMIN
xlo   <- rangex[1]
xhi   <- rangex[2]
width <- xhi - xlo
xfine <- seq(xlo, xhi, len=nfine)

#  evaluate target curve at fine mesh of values

y0fine <- eval.fd(xfine, y0fd)

#  set up indices of coefficients that will be modified in ACTIVE

wcoef1   <- wcoef[1,]
if (periodic) {
   active   <- 1:nbasis
   wcoef[1] <- 0
   shift    <- 0
} else {
   active <- 2:nbasis
}

#  initialize matrix Kmat defining penalty term

lambda <- WfdParobj$lambda
if (lambda > 0) {
   Lfdobj <- WfdParobj$Lfd
   Kmat <- getbasispenalty(wbasis, Lfdobj)
   ind  <- 2:ncoef
   Kmat <- lambda*Kmat[ind,ind]
} else {
   Kmat <- NULL
}

#  set up limits on coefficient sizes

climit <- 50*c(-rep(1,ncoef), rep(1,ncoef))

#  set up cell for storing basis function values

JMAX <- 15
basislist <- vector("list", JMAX)

yregcoef <- yfd$coefs

#  iterate through the curves

wcoefnew <- wcoef
if (dbglev == 0 && ncurve > 1) cat("Progress:  Each dot is a curve\n")

for (icurve in 1:ncurve) {
  if (dbglev == 0 && ncurve > 1) cat(".")
  if (dbglev >= 1 && ncurve > 1)
      cat(paste("\n\n-------  Curve ",icurve,"  --------\n"))
  if (ncurve == 1) {
    yfdi <- yfd
    Wfdi <- Wfd0
    cvec <- wcoef
  } else {
    Wfdi <- Wfd0[icurve]
    cvec <- wcoef[,icurve]
    if (nvar == 1) {
      yfdi <- yfd[icurve]
    } else {
      yfdi <- yfd[icurve,]
    }
  }

  #  evaluate curve to be registered at fine mesh

  yfine <- eval.fd(xfine, yfdi)

  #  evaluate objective function for starting coefficients

  #  first evaluate warping function and its derivative at fine mesh

  ffine  <-   monfn(xfine, Wfdi, basislist)
  Dffine <- mongrad(xfine, Wfdi, basislist)
  fmax   <- ffine[nfine]
  Dfmax  <- Dffine[nfine,]
  hfine  <- xlo + width*ffine/fmax
  Dhfine <- width*(fmax*Dffine - outer(ffine,Dfmax))/fmax^2
  hfine[1]     <- xlo
  hfine[nfine] <- xhi

  #  register curves given current Wfdi

  yregfdi <- regyfn(xfine, yfine, hfine, yfdi, Wfdi, periodic)

  #  compute initial criterion value and gradient

  Fstr <- regfngrad(xfine, y0fine, Dhfine, yregfdi, Wfdi, Kmat, periodic, crit)

  #  compute the initial expected Hessian

  if (crit == 2) {
     D2hwrtc <- monhess(xfine, Wfdi, basislist)
     D2fmax  <- D2hwrtc[nfine,]
     fmax2 <- fmax*fmax
     fmax3 <- fmax*fmax2
     m <- 1
     for (j in 2:nbasis) {
        m <- m + 1
        for (k in 2:j) {
           m <- m + 1
           D2hwrtc[,m] <- width*(2*ffine*Dfmax[j]*Dfmax[k]
                - fmax*(Dffine[,j]*Dfmax[k] + Dffine[,k]*Dfmax[j])
                + fmax2*D2hwrtc[,m] - ffine*fmax*D2fmax[m])/fmax3
        }
     }
  } else {
     D2hwrtc <- NULL
  }

  hessmat <- reghess(xfine, y0fine, Dhfine, D2hwrtc, yregfdi, Kmat, periodic, crit)

  #  evaluate the initial update vector for correcting the initial cvec

  result   <- linesearch(Fstr, hessmat, dbglev)
  deltac   <- result[[1]]
  cosangle <- result[[2]]
  #  initialize iteration status arrays

  iternum <- 0
  status <- c(iternum, Fstr$f, Fstr$norm)
  if (dbglev >= 1) {
        cat("\nIter.    Criterion   Grad Length")
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2],4))
        cat("      ")
        cat(round(status[3],4))
  }
  iterhist <- matrix(0,iterlim+1,length(status))
  iterhist[1,]  <- status
  if (iterlim == 0) break

  #  -------  Begin main iterations  -----------

  MAXSTEPITER <- 5
  MAXSTEP <- 100
  trial   <- 1
  reset   <- 0
  linemat <- matrix(0,3,5)
  cvecold <- cvec
  Foldstr <- Fstr
  dbgwrd  <- dbglev >= 2
  #  ---------------  beginning of optimization loop  -----------
  for (iter in 1:iterlim) {
      iternum <- iternum + 1
      #  set logical parameters
      dblwrd <- c(FALSE,FALSE)
      limwrd <- c(FALSE,FALSE)
      ind <- 0
      ips <- 0
      #  compute slope
      linemat[2,1] <- sum(deltac*Foldstr$grad)
      #  normalize search direction vector
      sdg          <- sqrt(sum(deltac^2))
      deltac       <- deltac/sdg
      linemat[2,1] <- linemat[2,1]/sdg
      # initialize line search vectors
      linemat[,1:4] <- matrix(c(0, linemat[2,1], Fstr$f),3,1) %*% matrix(1,1,4)
      stepiter  <- 0
      if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,1]),4)))
      }
      #  return with stop condition if initial slope is nonnegative
      if (linemat[2,1] >= 0) {
        if (dbglev >= 2) cat("\nInitial slope nonnegative.")
        ind <- 3
        break
      }
      #  return successfully if initial slope is very small
      if (linemat[2,1] >= -min(c(1e-3,conv))) {
        if (dbglev >= 2) cat("\nInitial slope too small")
        ind <- 0
        break
      }
      #  first step set to trial
      linemat[1,5]  <- trial
      #  ------------  begin line search iteration loop  ----------
      cvecnew <- cvec
      Wfdnewi <- Wfdi
      for (stepiter in 1:MAXSTEPITER) {
        #  check the step size and modify if limits exceeded
        result <- stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                   climit, active, dbgwrd)
        linemat[1,5] <- result[[1]]
        ind          <- result[[2]]
        limwrd       <- result[[3]]
        if (ind == 1) break    # break of limit hit twice in a row
        if (linemat[1,5] <= 1e-7) {
           #  Current step size too small  terminate
           if (dbglev >= 2)
             cat("\nStepsize too small: ", round(linemat[1,5],4))
           break
        }
        #  update parameter vector
        cvecnew <- cvec + linemat[1,5]*deltac
        #  compute new function value and gradient
        Wfdnewi[[1]] <- cvecnew
        #  first evaluate warping function and its derivative at fine mesh
        cvectmp <- cvecnew
        cvectmp[1] <- 0
        Wfdtmpi <- Wfdnewi
        Wfdtmpi[[1]] <- cvectmp
        ffine  <-    monfn(xfine, Wfdtmpi, basislist)
        Dffine <- mongrad(xfine, Wfdtmpi, basislist)
        fmax   <- ffine[nfine]
        Dfmax  <- Dffine[nfine,]
        hfine  <- xlo + width*ffine/fmax
        Dhfine <- width*(fmax*Dffine - outer(ffine,Dfmax))/fmax^2
        hfine[1]     <- xlo
        hfine[nfine] <- xhi
        #  register curves given current Wfdi
        yregfdi <- regyfn(xfine, yfine, hfine, yfdi, Wfdnewi, periodic)
        Fstr    <- regfngrad(xfine, y0fine, Dhfine, yregfdi, Wfdnewi, Kmat, periodic, crit)
        linemat[3,5] <- Fstr$f
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*Fstr$grad)
        if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,5]),4)))
        }
        #  compute next line search step, also testing for convergence
        result  <- stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbgwrd)
        linemat <- result[[1]]
        ips     <- result[[2]]
        ind     <- result[[3]]
        dblwrd  <- result[[4]]
        trial   <- linemat[1,5]
        #  ind == 0 implies convergence
        if (ind == 0 || ind == 5) break
     }
     #  ------------  end line search iteration loop  ----------
     cvec   <- cvecnew
     Wfdi   <- Wfdnewi
     #  test for function value made worse
     if (Fstr$f > Foldstr$f) {
        #  Function value worse  warn and terminate
        ier <- 1
        if (dbglev >= 2) {
          cat("Criterion increased, terminating iterations.\n")
          cat(paste("\n",round(c(Foldstr$f, Fstr$f),4)))
        }
        #  reset parameters and fit
        cvec   <- cvecold
        Wfdi[[1]] <- cvecold
        Fstr   <- Foldstr
        deltac <- -Fstr$grad
        if (dbglev > 2) {
          for (i in 1:nbasis) cat(cvec[i])
          cat("\n")
        }
        if (reset == 1) {
           #  This is the second time in a row that this
           #     has happened   quit
           if (dbglev >= 2) cat("Reset twice, terminating.\n")
            break
        } else {
           reset <- 1
        }
     } else {
        #  function value has not increased,  check for convergence
        if (abs(Foldstr$f-Fstr$f) < conv) {
           wcoef[,icurve]    <- cvec
           status <- c(iternum, Fstr$f, Fstr$norm)
           iterhist[iter+1,] <- status
           if (dbglev >= 1) {
              cat("\n")
              cat(iternum)
              cat("        ")
              cat(round(status[2],4))
              cat("      ")
              cat(round(status[3],4))
	        }
           break
        }
        #  update old parameter vectors and fit structure
        cvecold <- cvec
        Foldstr <- Fstr
        #  update the expected Hessian
        if (crit == 2) {
           cvectmp <- cvec
           cvectmp[1] <- 0
           Wfdtmpi[[1]] <- cvectmp
           D2hwrtc <- monhess(xfine, Wfdtmpi, basislist)
           D2fmax  <- D2hwrtc[nfine,]
           #  normalize 2nd derivative
           fmax2 <- fmax*fmax
           fmax3 <- fmax*fmax2
           m <- 1
           for (j in 2:nbasis) {
              m <- m + 1
              for (k in 2:j) {
                 m <- m + 1
                 D2hwrtc[,m] <- width*(2*ffine*Dfmax[j]*Dfmax[k]
                - fmax*(Dffine[,j]*Dfmax[k] + Dffine[,k]*Dfmax[j])
                + fmax2*D2hwrtc[,m] - ffine*fmax*D2fmax[m])/fmax3
              }
           }
        } else {
           D2hwrtc <- NULL
        }
        hessmat <- reghess(xfine, y0fine, Dhfine, D2hwrtc, yregfdi, Kmat, periodic, crit)
        #  update the line search direction vector
        result   <- linesearch(Fstr, hessmat, dbglev)
        deltac   <- result[[1]]
        cosangle <- result[[2]]
        reset <- 0
     }
     status <- c(iternum, Fstr$f, Fstr$norm)
     iterhist[iter+1,] <- status
     if (dbglev >= 1) {
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2],4))
        cat("      ")
        cat(round(status[3],4))
     }
   }
  #  ---------------  end of optimization loop  -----------
  wcoef[,icurve] <- cvec
  if (nvar == 1) {
     yregcoef[,icurve]  <- yregfdi$coefs
  } else {
     yregcoef[,icurve,] <- yregfdi$coefs
  }
}

#  --------------------   end of variable loop  -----------

#  create functional data objects for the registered curves

regfdnames <- yfd$fdnames
regfdnames[[3]] <- paste("Registered ",regfdnames[[3]])
ybasis  <- yfd$basis
regfd   <- fd(yregcoef, ybasis, regfdnames)

#  create functional data objects for the warping functions

if (periodic) {
  shift <- c(wcoef[1,])
  wcoef[1,] <- wcoef1
} else {
  shift <- rep(0,ncurve)
}
Wfd <- fd(wcoef, wbasis)

regstr <- list("regfd"=regfd, "Wfd"=Wfd, "shift"=shift)

return(regstr)
}

#  ----------------------------------------------------------------

regfngrad <- function(xfine, y0fine, Dhwrtc, yregfd, Wfd, Kmat, periodic, crit)
{
	#cat("\nregfngrad")
  y0dim <- dim(y0fine)
  if (length(y0dim) == 3) nvar <- y0dim[3] else nvar <- 1
  nfine <- length(xfine)
  cvec  <- Wfd$coefs
  ncvec <- length(cvec)
  onecoef <- matrix(1,1,ncvec)

  if (periodic) {
     Dhwrtc[,1] <- 1
  } else {
     Dhwrtc[,1] <- 0
  }
  yregmat  <- eval.fd(xfine, yregfd)
  Dyregmat <- eval.fd(xfine, yregfd, 1)
  #if (nvar > 1) {
  #	 y0fine   <- y0fine[,1,]
  #	 yregmat  <- yregmat[,1,]
  #	 Dyregmat <- Dyregmat[,1,]
  #}

  #  loop through variables computing function and gradient values

  Fval <- 0
  gvec <- matrix(0,ncvec,1)
  for (ivar in 1:nvar) {
    y0ivar  <-   y0fine[,ivar]
    ywrthi  <-  yregmat[,ivar]
    Dywrthi <- Dyregmat[,ivar]
    aa      <- mean(y0ivar^2)
    bb      <- mean(y0ivar*ywrthi)
    cc      <- mean(ywrthi^2)
    Dywrtc  <- (Dywrthi %*% onecoef)*Dhwrtc
    if (crit == 1) {
      res  <- y0ivar - ywrthi
      Fval <- Fval + aa - 2*bb + cc
      gvec <- gvec - 2*crossprod(Dywrtc, res)/nfine
    } else {
      ee   <- aa + cc
      ff   <- aa - cc
      dd   <- sqrt(ff^2 + 4*bb^2)
      Fval <- Fval + ee - dd
      Dbb  <- crossprod(Dywrtc, y0ivar)/nfine
      Dcc  <- 2.0 * crossprod(Dywrtc, ywrthi)/nfine
      Ddd  <- (4*bb*Dbb - ff*Dcc)/dd
      gvec <- gvec + (Dcc - Ddd)
    }
  }
  if (!is.null(Kmat)) {
     ind   <- 2:ncvec
     ctemp <- cvec[ind,1]
     Kctmp <- Kmat%*%ctemp
     Fval  <- Fval + t(ctemp)%*%Kctmp
     gvec[ind] <- gvec[ind] + 2*Kctmp
  }

#  set up FALSE structure containing function value and gradient
  Fstr <- list(f=0, grad=rep(0,ncvec), norm=0)
  Fstr$f    <- Fval
  Fstr$grad <- gvec
  #  do not modify initial coefficient for B-spline and Fourier bases
  if (!periodic)  Fstr$grad[1] <- 0
  Fstr$norm <- sqrt(sum(Fstr$grad^2))
  return(Fstr)
}

#  ---------------------------------------------------------------

reghess <- function(xfine, y0fine, Dhwrtc, D2hwrtc, yregfd, Kmat,
                           periodic, crit)
{
	#cat("\nreghess")
  y0dim <- dim(y0fine)
  if (length(y0dim) == 3) nvar <- y0dim[3] else nvar <- 1
  nfine   <- length(xfine)
  ncoef   <- dim(Dhwrtc)[2]
  onecoef <- matrix(1,1,ncoef)
  npair   <- ncoef*(ncoef+1)/2

  if (periodic) {
     Dhwrtc[,1] <- 1
  } else {
     Dhwrtc[,1] <- 0
  }
  yregmat  <- eval.fd(yregfd, xfine)
  Dyregmat <- eval.fd(yregfd, xfine, 1)
  if (nvar > 1) {
	 y0fine   <- y0fine[,1,]
	 yregmat  <- yregmat[,1,]
	 Dyregmat <- Dyregmat[,1,]
  }

  if (crit == 2) {
     D2yregmat <- eval.fd(yregfd, xfine, 2)
     if (nvar > 1) D2yregmat <- D2yregmat[,1,]
     if (periodic) {
        D2hwrtc[,1] <- 0
        for (j in 2:ncoef) {
           m <- j*(j-1)/2 + 1
           D2hwrtc[,m] <- Dhwrtc[,j]
        }
     } else {
        D2hwrtc[,1] <- 1
        for (j in 2:ncoef) {
           m <- j*(j-1)/2 + 1
           D2hwrtc[,m] <- 0
        }
     }
  }

  hessvec <- matrix(0,npair,1)
  for (ivar in 1:nvar) {
    y0i        <-   y0fine[,ivar]
    yregmati   <-  yregmat[,ivar]
    Dyregmati  <- Dyregmat[,ivar]
    Dywrtc <- ((Dyregmati %*% onecoef)*Dhwrtc)
    if (crit == 1) {
      hessmat <-  2*crossprod(Dywrtc, Dywrtc)/nfine
      m <- 0
       for (j in 1:ncoef) {
        for (k in 1:j) {
          m <- m + 1
          hessvec[m] <- hessvec[m] + hessmat[j,k]
        }
      }
    } else {
      D2yregmati <- D2yregmat[,ivar]
      aa     <- mean(y0i^2)
      bb     <- mean(y0i*yregmati)
      cc     <- mean(    yregmati^2)
      Dbb    <- crossprod(Dywrtc, y0i)/nfine
      Dcc    <- 2.0 * crossprod(Dywrtc, yregmati)/nfine
      D2bb   <- matrix(0,npair,1)
      D2cc   <- matrix(0,npair,1)
      crossprodmat <- matrix(0,nfine,npair)
      DyD2hmat     <- matrix(0,nfine,npair)
      m <- 0
      for (j in 1:ncoef) {
        for (k in 1:j) {
          m <- m + 1
          crossprodmat[,m] <- Dhwrtc[,j]*Dhwrtc[,k]*D2yregmati
          DyD2hmat[,m] <- Dyregmati*D2hwrtc[,m]
          temp <- crossprodmat[,m] + DyD2hmat[,m]
          D2bb[m] <- mean(y0i*temp)
          D2cc[m] <- 2*mean(yregmati*temp +
                     Dyregmati^2*Dhwrtc[,j]*Dhwrtc[,k])
        }
      }
      ee     <- aa + cc
      ff     <- aa - cc
      ffsq   <- ff*ff
      dd     <- sqrt(ffsq + 4*bb*bb)
      ddsq   <- dd*dd
      ddcu   <- ddsq*dd
      m <- 0
      for (j in 1:ncoef) {
        for (k in 1:j) {
          m <- m + 1
          hessvec[m] <- hessvec[m] + D2cc[m] -
            (4*Dbb[j]*Dbb[k] + 4*bb*D2bb[m] + Dcc[j]*Dcc[k] -
                   ff* D2cc[m])/dd +
            (4*bb*Dbb[j] - ff*Dcc[j])*(4*bb*Dbb[k] - ff*Dcc[k])/ddcu
        }
      }
    }
  }
  hessmat <- matrix(0,ncoef,ncoef)
  m <- 0
  for (j in 1:ncoef) {
    for (k in 1:j) {
      m <- m + 1
      hessmat[j,k] <- hessvec[m]
      hessmat[k,j] <- hessvec[m]
    }
  }
  if (!is.null(Kmat)) {
     ind <- 2:ncoef
     hessmat[ind,ind] <- hessmat[ind,ind] + 2*Kmat
  }
  if (!periodic) {
     hessmat[1,]  <- 0
     hessmat[,1]  <- 0
     hessmat[1,1] <- 1
  }
  return(hessmat)
}

#  ----------------------------------------------------------------

regyfn <- function(xfine, yfine, hfine, yfd, Wfd, periodic)
{
	#cat("\nregyfn")
coef  <- Wfd$coefs
shift <- coef[1]
coef[1] <- 0
Wfd[[1]] <- coef

if (all(coef == 0)) {
   if (periodic) {
      if (shift == 0) {
         yregfd <- yfd
         return(yregfd)
      }
   } else {
      yregfd <- yfd
      return(yregfd)
   }
}

#  Estimate inverse of warping function at fine mesh of values
#  28 dec 000
#  It makes no real difference which
#     interpolation method is used here.
#  Linear is faster and sure to be monotone.
#  Using WARPSMTH added nothing useful, and was abandoned.
nfine       <- length(xfine)
hinv        <- approx(hfine, xfine, xfine)$y
hinv[1]     <- xfine[1]
hinv[nfine] <- xfine[nfine]

#  carry out shift if period and shift != 0
basis  <- yfd$basis
rangex <- basis$rangeval
ydim <- dim(yfine)
#if (length(ydim) == 3) yfine <- yfine[,1,]
if (periodic & shift != 0) yfine <- shifty(xfine, yfine, shift)
#  make FD object out of Y
ycoef  <- project.basis(yfine, hinv, basis, 1)
yregfd <- fd(ycoef, basis)
return(yregfd)
}

#  ----------------------------------------------------------------

linesearch <- function(Fstr, hessmat, dbglev)
{
deltac   <- -solve(hessmat,Fstr$grad)
cosangle <- -sum(Fstr$grad*deltac)/sqrt(sum(Fstr$grad^2)*sum(deltac^2))
if (dbglev >= 2) cat(paste("\nCos(angle) = ",round(cosangle,2)))
if (cosangle < 1e-7) {
   if (dbglev >=2) cat("\nangle negative")
   deltac <- -Fstr$grad
}
return(list(deltac, cosangle))
}

#  ---------------------------------------------------------------------

shifty <- function(x, y, shift)
{
#SHIFTY estimates value of Y for periodic data for
#       X shifted by amount SHIFT.
#  It is assumed that X spans interval over which functionis periodic.
#  Last modified 6 February 2001

ydim <- dim(y)
if (is.null(ydim)) ydim <- 1
if (length(ydim) > 3) stop("Y has more than three dimensions")

if (shift == 0) {
   yshift <- y
   return(yshift)
}

n   <- ydim[1]
xlo <- min(x)
xhi <- max(x)
wid <- xhi - xlo
if (shift > 0) {
   while (shift > xhi)  shift <- shift - wid
   ind <- 2:n
   x2  <- c(x, x[ind]+wid)
   xshift <- x + shift
   if (length(ydim) == 1) {
	  y2 <- c(y, y[ind])
      yshift <- approx(x2, y2, xshift)$y
   }
   if (length(ydim) == 2) {
	   nvar <- ydim[2]
	   yshift <- matrix(0,n,nvar)
      for (ivar in 1:nvar) {
         y2 <- c(y[,ivar], y[ind,ivar])
         yshift[,ivar] <- approx(x2, y2, xshift)$y
      }
   }
   if (length(ydim) == 3) {
	   nrep <- ydim[2]
	   nvar <- ydim[3]
      yshift <- array(0,c(n,nrep,nvar))
      for (irep in 1:nrep) for (ivar in 1:nvar) {
         y2 <- c(y[,irep,ivar], y[ind,irep,ivar])
         yshift[,irep,ivar] <- approx(x2, y2, xshift)$y
      }
   }
} else {
   while (shift < xlo - wid) shift <- shift + wid
   ind <- 1:(n-1)
   x2 <- c(x[ind]-wid, x)
   xshift <- x + shift
   if (length(ydim) == 1) {
      y2 <- c(y[ind], y)
      yshift <- approx(x2, y2, xshift)$y
   }
   if (length(ydim) == 2) {
	   nvar <- ydim[2]
	   yshift <- matrix(0,n,nvar)
	   for (ivar in 1:nvar) {
		   y2 <- c(y[ind,ivar],y[,ivar])
		   yshift[,ivar] <- approx(x2, y2, xshift)$y
	   }
   }
   if (length(ydim) == 3) {
	   nrep <- ydim[2]
	   nvar <- ydim[3]
      yshift <- array(0, c(n,nrep,nvar))
      for (irep in 1:nrep) for (ivar in 1:nvar) {
         y2 <- c(y[ind,irep,ivar], y[,irep,ivar])
         yshift[,irep,ivar] <- approx(x2, y2, xshift)$y
      }
   }
}
return(yshift)
}

