#  ----------------------------------------------------------------------------

eval.fd <- function(evalarg, fdobj, Lfdobj=0) {

#  EVAL_FD evaluates a functional data observation at argument
#  values EVALARG.
#
#  LFDOBJ is a functional data object defining the order m
#  HOMOGENEOUS linear differential operator of the form
#  Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) +
#          \exp[w_m(t)] D^m x(t)
#
#  Arguments:
#  EVALARG ... A vector of values at which all functions are to
#              evaluated.
#  FDOBJ   ... Functional data object
#  LFDOBJ  ... A linear differential operator object
#              applied to the functions before they are evaluated.

#  Note that the first two arguments may be interchanged.

#  Returns:  An array of function values corresponding to the evaluation
#              arguments in EVALARG

#  Last modified 3 January 2008 by Jim Ramsay

#  Check LFDOBJ

Lfdobj <- int2Lfd(Lfdobj)

#  Exchange the first two arguments if the first is an FD object
#    and the second numeric

if (inherits(fdobj, "numeric") && inherits(evalarg, "fd")) {
  temp    <- fdobj
  fdobj   <- evalarg
  evalarg <- temp
}

#  check EVALARG

if (!(is.numeric(evalarg))) stop("Argument EVALARG is not numeric.")

evaldim <- dim(evalarg)
if (!(length(evaldim) < 3)) stop(
   "Argument EVALARG is not a vector or a matrix.")

#  check FDOBJ

if (!(inherits(fdobj, "fd"))) stop(
     "Argument FD is not a functional data object.")

#  Extract information about the basis

basisobj <- fdobj$basis
nbasis   <- basisobj$nbasis
rangeval <- basisobj$rangeval
onerow   <- rep(1,nbasis)

temp <- c(evalarg)
temp <- temp[!(is.na(temp))]
EPS  <- 1e-14
if (min(temp) < rangeval[1]-EPS || max(temp) > rangeval[2]+EPS) {
  warning(paste(
    "Values in argument EVALARG are outside of permitted range,",
    "and will be ignored."))
    print(c(rangeval[1]-min(temp), max(temp) - rangeval[2]))
}

#  get maximum number of evaluation values

if (is.vector(evalarg)) {
  n <- length(evalarg)
} else {
  n <- evaldim[1]
}

#  Set up coefficient array for FD

coef  <- fdobj$coefs
coefd <- dim(coef)
ndim  <- length(coefd)
if (ndim <= 1) nrep <- 1 else nrep <- coefd[2]
if (ndim <= 2) nvar <- 1 else nvar <- coefd[3]

#  Set up array for function values

if (ndim <= 2) evalarray <- matrix(0,n,nrep)
else           evalarray <- array(0,c(n,nrep,nvar))
if (ndim == 2) dimnames(evalarray) <- list(NULL,dimnames(coef)[[2]])
if (ndim == 3) dimnames(evalarray) <- list(NULL,dimnames(coef)[[2]],
                                                dimnames(coef)[[3]])

#  Case where EVALARG is a vector of values to be used for all curves

if (is.vector(evalarg)) {

    evalarg[evalarg < rangeval[1]-1e-10] <- NA
    evalarg[evalarg > rangeval[2]+1e-10] <- NA
    basismat <- eval.basis(evalarg, basisobj, Lfdobj)

    #  evaluate the functions at arguments in EVALARG

    if (ndim <= 2) {
      evalarray <- basismat %*% coef
    } else {
       evalarray <- array(0,c(n,nrep,nvar))
       for (ivar in 1:nvar) evalarray[,,ivar] <- basismat %*% coef[,,ivar]
    }

} else {

  #  case of evaluation values varying from curve to curve

  for (i in 1:nrep) {
    evalargi <- evalarg[,i]
       if (all(is.na(evalargi))) stop(
            paste("All values are NA for replication",i))

    index    <- !(is.na(evalargi) | evalargi < rangeval[1] |
                                    evalargi > rangeval[2])
    evalargi <- evalargi[index]
       basismat <- eval.basis(evalargi, basisobj, Lfdobj)

       #  evaluate the functions at arguments in EVALARG

       if (ndim == 2) {
           evalarray[  index, i] <- basismat %*% coef[,i]
           evalarray[!(index),i] <- NA
       }
       if (ndim == 3) {
           for (ivar in 1:nvar) {
             evalarray[   index,i,ivar] <- basismat %*% coef[,i,ivar]
               evalarray[!(index),i,ivar] <- NA
           }
       }
  }

}

return(evalarray)

}

