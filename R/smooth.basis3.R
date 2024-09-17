smooth.basis3 <- function(argvals=array(1:n,c(n,N,M)), y, fdParobj,
                          wtvec=NULL,   fdnames=NULL, covariates=NULL,
                          method="chol", dfscale=1, returnMatrix=FALSE)
{
  ##
  ## 1.  check dimensions of argval and y
  ##
  
  dimy <- dim(y)
  ndy <- length(dimy)
  n   <- dimy[1]
  N   <- dimy[2]
  M   <- dimy[3]
  if (ndy < 3)stop("length(dim(y)) must be 3  is ", ndy)
  if (any(dima != dimy)) {
    stop("dim(argvals) = ", paste(dima, collapse=", "),
         " != dim(y) = ", paste(dimy, collapse=", "))
  }
  
  dima <- dim(argvals)
  nda  <- length(dima)
  if (nda < 3) stop("length(dim(argvals)) must be 3  is ", nda)
  
  ##
  ## 2.  Call smooth.basis2 repeatedly
  ##
  #  2.1.  argvals[, , 1]
  sb1 <- smooth.basis2(argvals[, , 1], y=y[, , 1], fdParobj=fdParobj,
                       wtvec=wtvec,   fdnames=fdnames,
                       covariates=covariates,
                       method=method, dfscale=dfscale,
                       returnMatrix=returnMatrix)
  #  2.2.  set up output object
  coef1 <- sb1$fd$coefs
  dimc1 <- dim(coef1)
  dimc  <- c(dimc1[1], dimy[-1])
  coefs <- array(NA, dim=dimc)
  argNames <- dimnames(argvals)
  yNames   <- dimnames(y)
  c1Names  <- dimnames(coef1)
  cNames   <- vector("list", 3)
  if (!is.null(c1Names[[1]])) cNames[[1]] <- c1Names[[1]]
  if (!is.null(yNames[[2]]))  cNames[[2]] <- yNames[[2]]
  if (!is.null(yNames[[3]]))  cNames[[3]] <- yNames[[3]]
  dimnames(coefs) <- cNames
  if (!is.null(covariates)) {
    q <- dim(covariates)[2]
    beta. <- array(0,c(q,dimy[2],dimy[3]))
    beta.[,,1] <- sb1$beta
  } else {
    beta. <- NULL
  }
  #
  for (i in seq(2, length=dimy[3]-1)) {
    sbi <- smooth.basis2(argvals[,,i], y=y[,,i], fdParobj=fdParobj,
                         wtvec=wtvec,   fdnames=fdnames,
                         covariates=covariates,
                         method=method, dfscale=dfscale,
                         returnMatrix=returnMatrix)
    coefs[,,i] <- sbi$fd$coefs
    if (!is.null(covariates)) {
      beta.[,,i] <- sbi$beta
    }
  }
  if (is.null(fdnames)) {
    fdnames <- list(time=NULL, reps=NULL, values=NULL)
    if (!is.null(yNames[[1]])) {
      fdnames[[1]] <- yNames[[1]]
    } else {
      if (!is.null(argNames[[1]]))
        fdnames[[1]] <- argNames[[1]]
    }
    if (!is.null(yNames[[2]])) {
      fdnames[[2]] <- yNames[[2]]
    } else {
      if (!is.null(argNames[[2]]))
        fdnames[[2]] <- argNames[[2]]
    }
    if (!is.null(yNames[[3]])) {
      fdnames[[3]] <- yNames[[3]]
    } else {
      if (!is.null(argNames[[3]]))
        fdnames[[3]] <- argNames[[3]]
    }
  }
  ##
  ## 3.  done
  ##
  sb <- sb1
  sb$fd$coefs   <- coefs
  sb$fd$fdnames <- fdnames
  sb$beta       <- beta.
  sb
  
}

