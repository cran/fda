smooth.basis2 <- function(argvals=matrix(1:n,n,N), y, fdParobj,
                          wtvec=NULL,   fdnames=NULL, covariates=NULL,
                          method="chol", dfscale=1, returnMatrix=FALSE) {
  ##
  ## 1.  number of  dimensions of y = 2 or 3?
  ##
  dimy     <- dim(y)
  ndy      <- length(dimy)
  n        <- dimy[1]
  N        <- dimy[2]
  ynames   <- dimnames(y)
  argNames <- dimnames(argvals)
  ##
  ## 2.  ndy == 2
  ##
  if (ndy < 3) {
    #  2.1.  Start by smoothing first record using argvals[, 1]
    
    sb1 <- smooth.basis1(argvals[, 1], y=y[, 1], fdParobj=fdParobj,
                         wtvec=wtvec,   fdnames=fdnames,
                         covariates=covariates,
                         method=method, dfscale=dfscale,
                         returnMatrix=returnMatrix)
    
    #  2.2.  set up output object
    dimc1   <- dim(sb1$fd$coefs)
    dimc    <- c(dimc1[1], dimy[-1])
    coefs   <- array(NA, dim=dimc)
    c1names <- dimnames(sb1$fd$coefs)
    cNames  <- vector("list", 2)
    if (!is.null(c1names[[1]])) cNames[[1]] <- c1names[[1]]
    if (!is.null(ynames[[2]]))  cNames[[2]] <- ynames[[2]]
    dimnames(coefs) <- cNames
    coefs[, 1] <- sb1$fd$coefs
    if (!is.null(covariates)) {
      q <- dim(covariates)[2]
      beta. <- matrix(0,q,dimy[2])
      beta.[,1] <- sb1$beta
    } else {
      beta. <- NULL
    }
    #   now loop through remaining records, smoothing each in term,
    #   using argvals[,1]
    for (i in seq(2, length=dimy[2]-1)) {
      
      sbi <- smooth.basis1(argvals[, i], y=y[, i], fdParobj=fdParobj,
                           wtvec=wtvec,   fdnames=fdnames,
                           covariates=covariates,
                           method=method, dfscale=dfscale,
                           returnMatrix=returnMatrix)
      
      coefs[, i] <- sbi$fd$coefs
      if (!is.null(covariates)) {
        beta.[,i] <- sbi$beta
      }
    }
    if (is.null(fdnames)) {
      fdnames <- sb1$fdnames
      if (is.null(fdnames))
        fdnames <- list(time=NULL, reps=NULL, values="value")
      valueChk <- ((length(fdnames$values)==1)
                   && (fdnames$values=="value")
                   && (length(fdnames$reps)==1)
                   && (!is.null(ynames[[2]])) )
      if (valueChk)fdnames$values <- fdnames$reps
      if (!is.null(ynames[[2]]))
        fdnames[[2]] <- ynames[[2]]
    }
  } else {
    ##
    ## 3.  ndy == 3
    ##
    #  3.1.  argvals[, 1]
    sb1 <- smooth.basis1(argvals[, 1], y=y[, 1, ], fdParobj=fdParobj,
                         wtvec=wtvec,   fdnames=fdnames,
                         covariates=covariates,
                         method=method, dfscale=dfscale,
                         returnMatrix=returnMatrix)
    #  3.2.  set up output object
    coef1 <- sb1$fd$coefs
    dimc1 <- dim(coef1)
    dimc <- c(dimc1[1], dimy[-1])
    coefs <- array(NA, dim=dimc)
    yNames <- dimnames(y)
    c1Names <- dimnames(coef1)
    cNames <- vector("list", 3)
    if (!is.null(c1Names[[1]]))  cNames[[1]] <- c1Names[[1]]
    if (!is.null(yNames[[2]]))   cNames[[2]] <- yNames[[2]]
    if (is.null(c1Names[[2]])) {
      if (!is.null(yNames[[3]])) cNames[[3]] <- yNames[[3]]
    } else {
      cNames[[3]] <- c1Names[[2]]
    }
    dimnames(coefs) <- cNames
    coefs[, 1, ] <- coef1
    if (!is.null(covariates)) {
      q <- dim(covariates)[2]
      beta. <- array(0,c(q,dimy[2],dimy[3]))
      beta.[,,1] <- sb1$beta
    } else {
      beta. <- NULL
    }
    #
    for (i in seq(2, length=dimy[2]-1)) {
      sbi <- smooth.basis1(argvals[, i], y=y[, i, ], fdParobj=fdParobj,
                           wtvec=wtvec,   fdnames=fdnames,
                           covariates=covariates,
                           method=method, dfscale=dfscale)
      coefs[, i, ] <- sbi$fd$coefs
      if (!is.null(covariates)) {
        beta.[,,i] <- sbi$beta
      } else {
        beta. <- NULL
      }
    }
    if (is.null(fdnames)) {
      fdnames <- sb1$fdnames
      if (is.null(fdnames)) {
        fdnames <- list(time=NULL, reps=NULL, values=NULL)
        if (!is.null(argNames[[1]])) {
          fdnames[[1]] <- argNames[[1]]
        } else {
          fdnames[[1]] <- ynames[[1]]
        }
        if (!is.null(ynames[[2]]))fdnames[[2]] <- ynames[[2]]
        if (!is.null(ynames[[3]]))fdnames[[3]] <- ynames[[3]]
      }
    }
  }
  ##
  ## 4.  done
  ##
  sb <- sb1
  sb$beta       <- beta.
  sb$fd$coefs   <- coefs
  sb$fd$fdnames <- fdnames
  sb
  
}
