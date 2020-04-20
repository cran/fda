glm.fda <- function(basismat, y, family, lamRmat, wtvec=NULL, 
                    bvec0=NULL, addterm=NULL) {
  #GLM.FDA Fits a generalized linear model with regularization. 
  #  This function is called by function smooth.GLM
  #  Arguments
  #
  #  BASISMAT An N by NBASIS matrix of values of basis functions 
  #  Y        May be
  #                a N by NCURVE matrix of data to be fitted
  #                or, in the binomial case with local sample sizes M.i,
  #                a list array of length 2, the first of which cantains
  #                the matrix above containing observed frequencies,
  #                and the second of which contains the corresponding 
  #                sample sizes.  Note that in the binary or Bernoulli case,
  #                Y may be a matrix of 1"s and 0"s and the M"s are
  #                taken to be 1"s.
  #  FAMILY   A string indicating which of the five GLM family members
  #              is assumed
  #              "normal" or "gaussian" or "Gaussian"
  #              "binomial" or "binary" or "Bernoulli"
  #              "poisson"
  #              "gamma"
  #              "inverse gaussian" or "inverse Gaussian"
  #              or a list array of length(N) with each list containing
  #              a specification of the GLM family of a single observation.
  #  LAMRMAT  a \lambda R, that is, a roughness penalty matrix R of 
  #              order equal to the number of basis functions used or number
  #              of columns of basismat multiplied by a scalar roughness 
  #              penalty parameter \lambda
  #  wtvec    a vector of prior weights, such as the inverses of the
  #              relative variance of each observation.
  #  BVEC0    starting values for regresscion coefficients
  #  ADDTERM  a addterm with a coefficient fixed at 1.0.
  #
  #  Returns
  #  BVEC      Final estimate of coefficients
  #  DEVIANCE  Deviance values
  #
  #   Last modified 17 May 2018 by Jim Ramsay
  
  #--------------------------------------------------------------------------
  #                    Check arguments
  #--------------------------------------------------------------------------
  
  #  dimensions of basismat
  
  basismatDim <- dim(basismat)
  n       <- basismatDim[1]
  nbasis  <- basismatDim[2]
  if (is.list(y)) {
    yDim <- dim(as.matrix(y[[1]]))
    ntmp    <- yDim[1]
    ncurve  <- yDim[2]
  } else { 
    y    <- as.matrix(y)
    yDim <- dim(y)
    ntmp    <- yDim[1]
    ncurve  <- yDim[2]
  }
  if (n != ntmp) {
    stop("basismat and y do not have the same number of rows.")
  }
  
  #  define default weight vector wtvec and check for positivity
  
  if (is.null(wtvec)) {
    wtvec <- matrix(1,n,1)
  }
  
  if (any(wtvec <= 0)) {
    stop("Non-positive values of wtvec found.")
  }
  
  #--------------------------------------------------------------------------
  #  Process y and define anonymous functions according to the 
  #  distribution of y
  #     devFn    the deviance or loss function, 
  #                 called after convergence is achieved
  #     stdFn    the scale factor multiplying D eta
  #                 called second inside loop        
  #     linkFn   link function, eta <- linkFn(mu),
  #                 called prior to loop, maps data space into real line
  #     DlinkFn  derivative of the link function wrt to mu
  #                 called first inside loop
  #     IlinkFn  the inverse of the link function IlinkFn[eta] <- mu,
  #                 called last inside loop, maps eta into data space
  # Then set a starting value for the mean mu, avoiding boundary values.
  #--------------------------------------------------------------------------
  
  M <- NULL
  if (is.character(family)) {
    #  --------------------------------------------------------------------
    #    All observations are in the same family, family is a string
    #  --------------------------------------------------------------------
    if (!(family == "normal"   ||
          family == "binomial" ||
          family == "poisson"  ||
          family == "gamma"    ||
          family == "inverse gaussian")) {
      stop("The distribution is not valid.")
    }
    if (family == "normal") {
      #  Note  y can be any real number, no restrictions
      devFn   <- function(mu,y) (y - mu)^2
      stdFn   <- function(mu)  matrix(1,dim(mu))
      linkFn  <- function(mu)  mu
      DlinkFn <- function(mu)  matrix(1,dim(mu))
      IlinkFn <- function(eta) eta
      mu      <- y
    } 
    #  --------------------------------------------------------------------
    if (family == "binomial") {
      if (is.numeric(y)) {
        #  If y a matrix, M is taken to be 1 (set below)
        #  and it must be a binary matrix containing only 0"s and 1"s
        if (any(y < 0 | y > 1)) {
          stop(c("For binomial case, y a single column but ", 
                 " contains values other than 0 or 1."))
        }
        M <- matrix(1,n,ncurve)
      } else {
        if (is.list(y) && length(y) == 2) {
          #  If y is a list array of length 2, then first list 
          #  contains a matrix containing the number of successes and 
          #  the second list either contains a matrix of the same 
          #  size as the matrix in y{1} or a single positive 
          #  integer.  
          #  These values or this value is the number of trials M
          #  for a binomial or bernoulli distribution.
          #  M must be a positive integer.
          Freq <- y[[1]]
          M    <- y[[2]]
          if (length(M) == 1) {
            M <- M*matrix(1,n,ncurve)
          }
          if (!all(dim(M) == dim(Freq))) {
            stop(c("FAMILY is binomial and matrix M is not the same ", 
                   "size as matrix FREQ"))
          }
          if (any(M < 0)) {
            stop(c("FAMILY is binomial and one or more values in M ", 
                   "have nonpositive values"))
          }
          if (any(any(floor(M) != M))) {
            stop(c("FAMILY is binomial and one or more values in M ", 
                   "have noninteger values."))
          }
          #  Redefine y is the proportion of sucesses
          y <- Freq/M
        } else {
          stop(c("FAMILY is binomial and y has incorrect dimensions ", 
                 " or is of wrong type."))
        }
        devFn   <- function(mu,y) 2*M*(y*log((y+(y==0))/mu) + 
                                            (1-y)*log((1-y+(y==1))/(1-mu)))
        stdFn   <- function(mu)  sqrt(mu*(1-mu)/M)
        linkFn  <- function(mu)   log(mu/(1-mu))
        DlinkFn <- function(mu)    1/(mu*(1-mu))
        loBnd   <- -16
        upBnd   <- -loBnd
        IlinkFn <- function(eta) 1/(1 + exp(-constrain(eta,loBnd,upBnd)))
        mu      <- (M*y + 0.5)/(M + 1)
      }
    }
    #  --------------------------------------------------------------------
    if (family == "poisson") {
      #  Note y must not contain negative numbers
      if (any(y < 0)) {
        stop("FAMILY is poisson and y contains negative values")
      }
      devFn   <- function(mu,y) 2*(y*(log((y+(y==0))/mu)) - 
                                     (y - mu))
      stdFn   <- function(mu)  sqrt(mu)
      linkFn  <- function(mu)   log(mu)
      DlinkFn <- function(mu)     1/mu
      loBnd   <- -16
      upBnd   <- -loBnd
      IlinkFn <- function(eta) exp(constrain(eta,loBnd,upBnd))
      mu      <- y + 0.25
    }
    #  --------------------------------------------------------------------
    if (family == "gamma") {
      #  Note  y must contain only positive numbers
      if (any(y <= 0)) {
        stop("FAMILY is gamma and Y contains nonpositive values")
      }
      devFn   <- function(mu,y) 2*(-log(y/mu) + (y - mu)/mu)
      stdFn   <- function(mu)    mu
      linkFn  <- function(mu)  1/mu
      DlinkFn <- function(mu) -1/mu^2
      loBnd   <- -16
      upBnd   <- 1/loBnd
      IlinkFn <- function(eta) 1/constrain(eta,loBnd,upBnd)
      mu      <- max(y, eps)
    }
    #  --------------------------------------------------------------------
    if (family == "inverse gaussian") {
      #  Note  y must contain only positive numbers
      if (any(y <= 0)) {
        stop(c("FAMILY is inverse gaussian and Y contains ", 
               "nonpositive values"))
      }
      devFn   <- function(mu,y) ((y - mu)/mu)^2/ y
      stdFn   <- function(mu)  mu^(3/2)
      loBnd   <- -8
      upBnd   <- 1/loBnd
      linkFn  <- function(mu)  constrain(mu,loBnd,upBnd)^(-2)
      DlinkFn <- function(mu)  -2*mu^(-3)
      IlinkFn <- function(eta) constrain(eta,loBnd,upBnd)^(-1/2)
      mu      <- y
    }
  }
  # } else {
  #   if (is.list(family) && length(family) == n) {
  #     #  --------------------------------------------------------------------
  #     #    Observations can be in different families, family is a list array.
  #     #  --------------------------------------------------------------------
  #     mu      <- matrix(0,n,1)
  #     loBnd   <- matrix(0,n,1)
  #     upBnd   <- matrix(0,n,1)
  #     devFn   <- vector("list",n)
  #     stdFn   <- vector("list",n)
  #     linkFn  <- vector("list",n)
  #     DlinkFn <- vector("list",n)
  #     IlinkFn <- vector("list",n)
  #     #  Dealing with the presence of some binomial observations y has
  #     #  to be a list with n rows and 2 columns for all data.  Ugh!
  #     binomwrd <- is.list(y) && all(dim(y) == c(n,2))
  #   }
  #   for (i in 1:n) {
  #     familyi <- family[[i]]
  #     if (!is.character(familyi)) {
  #       stop("A distribution specification is not a string.")
  #     }
  #     if (family == "normal") {
  #       #  Note  y can be any real number, no restrictions
  #       devFn[[i]]   <- function(mu,y) (y - mu)^2
  #       stdFn[[i]]   <- function(mu)  matrix(1,dim(mu))
  #       linkFn[[i]]  <- function(mu)  mu
  #       DlinkFn[[i]] <- function(mu)  matrix(1,dim(mu))
  #       IlinkFn[[i]] <- function(eta) eta
  #       mu[i,] <- y[i,]
  #     }
  #     if (family == "binomial") {
  #       if (all(isnumeric(y[i,]))) {
  #         #  If y a matrix, M is taken to be 1 (set below)
  #         #  and it must be a binary matrix containing only 
  #         #0"s and 1"s
  #         if (any(y[i,] < 0 | y[i,] > 1)) {
  #           stop(c("For binomial case, y a single column but ", 
  #                  " contains values other than 0 or 1."))
  #         }
  #       } else {
  #         if (binomwrd) {
  #           Freqi <- y[[i,1]]
  #           Mi    <- y[[i,2]]
  #           if (length(Mi) == 1) {
  #             Mi <- Mi*matrix(1,1,ncurve)
  #           }
  #           if (!all(dim(Mi) == dim(Freqi))) {
  #             stop(paste("FAMILY is binomial and matrix M is not the same ", 
  #                        "dim as matrix FREQ"))
  #           }
  #           if (any(any(Mi < 0))) {
  #             stop(c("FAMILY is binomial and one or more values in M ", 
  #                    "have nonpositive values"))
  #           }
  #           if (any(any(floor(Mi) != Mi))) {
  #             stop(paste("FAMILY is binomial and one or more values in M ", 
  #                        "have noninteger values."))
  #           }
  #           #  Redefine y is the proportion of sucesses
  #           y[i,] <- (Freqi/Mi)
  #         } else {
  #           stop(paste("FAMILY is binomial and y has incorrect dimensions ", 
  #                      " or is of wrong type."))
  #         }
  #         devFn[[i]]   <- function(mu,y) 2*M*(y*log((y+(y==0))/mu) + 
  #                                               (1-y)*log((1-y+(y==1))/(1-mu)))
  #         stdFn[[i]]   <- function(mu)  sqrt(mu*(1-mu)/M)
  #         linkFn[[i]]  <- function(mu)   log(mu/(1-mu))
  #         DlinkFn[[i]] <- function(mu)    1/(mu*(1-mu))
  #         loBnd[i]   <- log(eps)
  #         upBnd[i]   <- -loBnd[i]
  #         IlinkFn[[i]] <- function(eta) 1/(1 + exp(-constrain(eta,loBnd,upBnd)))
  #         mu[i]      <- (M[i]*y[i] + 0.5)/(M[i] + 1)
  #       }
  #       if (family == "gamma") {
  #         #  Note  y must contain only positive numbers
  #         if (any(y[i] <= 0)) {
  #           stop("FAMILY is gamma and Y contains nonpositive values")
  #         }
  #         devFn[[i]]   <- function(mu,y) 2*(-log(y/mu) + (y - mu)/mu)
  #         stdFn[[i]]   <- function(mu)    mu
  #         linkFn[[i]]  <- function(mu)  1/mu
  #         DlinkFn[[i]] <- function(mu) -1/mu^2
  #         loBnd[i]   <- eps
  #         upBnd[i]   <- 1/loBnd[i]
  #         IlinkFn[[i]] <- function(eta) 1/constrain(eta,loBnd,upBnd)
  #         mu[i,]    <- max(y[i,], eps)
  #       }
  #       if (family == "inverse gaussian") {
  #         #  Note  y must contain only positive numbers
  #         if (any(y[i,] <= 0)) {
  #           stop(c("FAMILY is inverse gaussian and Y contains ", 
  #                  "nonpositive values"))
  #         }
  #         devFn[[i]]   <- function(mu,y) ((y - mu)/mu)^2/ y
  #         stdFn[[i]]   <- function(mu)  mu^(3/2)
  #         loBnd[i]   <- eps^(1/2)
  #         upBnd[i]   <- 1/loBnd[i]
  #         linkFn[[i]]  <- function(mu)  constrain(mu,loBnd,upBnd)^(-2)
  #         DlinkFn[[i]] <- function(mu)  -2*mu^(-3)
  #         IlinkFn[[i]] <- function(eta) constrain(eta,loBnd,upBnd)^(-1/2)
  #         mu[i,]    <- y[i,]
  #       }
  #     }
  #   }
  
  #--------------------------------------------------------------------------
  #                   Initialize mu and eta from y.
  #--------------------------------------------------------------------------
  
  # compute eta <- E(y) from mu
  
  if (is.character(family)) {
    eta <- linkFn(mu)
  # } else {
  #   eta  <- matrix(0,n,nurve)
  #   Deta <- matrix(0,n,nurve)
  #   stdm <- matrix(0,n,nurve)
  #   for (i in 1:n) {
  #     linkFni  <- linkFn[[i]]
  #     eta[i,] <- linkFni(mu[i,])
  #   }
  }
  
  #--------------------------------------------------------------------------
  #                        Set up for iterations
  #--------------------------------------------------------------------------
  
  iter     <- 0
  iterLim  <- 100
  seps     <- sqrt(eps)
  convcrit <- 1e-6
  sqrtwt   <- sqrt(wtvec)
  
  #  set up starting value bvec0 if required
  
  if (is.null(bvec0)) {
    bvec0 <- matrix(0,nbasis,ncurve)
  }
  bvec <- bvec0
  
  # Enforce limits on mu to guard against an inverse linkFn that doesn"t map 
  # into the support of the distribution.
  
  if (family == "binomial") {
    # mu is a probability, so order one is the natural scale, and eps is a
    # reasonable lower limit on that scale (plus it"s symmetric).
    eps <- 1e-16
    muLims <- c(eps, 1-eps)
  }
  if (family == "poisson" || family == "gamma" || family == "inverse gaussian") {
    # Here we don"t know the natural scale for mu, so make the lower limit
    # small.  This choice keeps mu^4 from underflowing.  No upper limit.
    muLims <- 1e-4
  }
  
  #--------------------------------------------------------------------------
  #                       Start of GLM iteration loop
  #--------------------------------------------------------------------------
  
  while (iter <= iterLim) {
    iter <- iter+1
    
    # Compute adjusted dep}ent variable for least squares fit
    
    if (is.character(family)) {
      Deta <- DlinkFn(mu)
      stdm <- stdFn(mu)
    # } else {
    #   for (i in 1:n) {
    #     DlinkFni  <- DlinkFn[[i]]
    #     stdFni    <- stdFn[[i]]
    #     mui       <- mu[i,]
    #     Deta[i,]  <- DlinkFni(mui)
    #     stdm[i,]  <- stdFni(mui)
    #   }
    }
    Zvec <- eta + (y - mu)*Deta
    
    # Compute IRLS weights the inverse of the variance function
    
    sqrtw <- (sqrtwt %*% matrix(1,1,ncurve))/(abs(Deta)*stdm)
    
    # Compute coefficient estimates for this iteration - the IRLS step
    
    bvec.old   <- bvec
    if (!is.null(addterm)) {
      ytmp <- Zvec - addterm
    } else {
      ytmp <- Zvec
    }
    yw   <- ytmp*sqrtw
    basismatw   <- basismat*(sqrtwt %*% matrix(1,1,nbasis))
    if (is.null(lamRmat)) {
      Mmat <- crossprod(basismatw)
    } else {
      Mmat <- crossprod(basismatw) + lamRmat
    }
    bvec    <- solve(Mmat,t(basismatw)) %*% yw
    if (!is.null(addterm)) {
      eta <- basismat %*% bvec + addterm
    } else {
      eta <- basismat %*% bvec
    }
    if (is.character(family)) {
      mu <- IlinkFn(eta)
    # } else {
    #   for (i in 1:n) {
    #     IlinkFni <- IlinkFn[[i]]
    #     mu[i,]   <- IlinkFni(eta[i,])
    #   }
    }
    
    # Force mean in bounds, in case the linkFn function is faulty
    
    if (is.character(family)) {
      if (family == "binomial") {
        if (any(mu < muLims[1] | muLims[2] < mu)) {
          for (j in 1:n) {
            mu[,j] <- max(min(mu[,j],muLims[2]),muLims[1])
          }
        }
      }
      if (family == "poisson" || 
          family == "gamma" || 
          family == "inverse gaussian") {
        if (any(mu < muLims[1])) {
          for (j in 1:n) {
            mu[j] <- max(mu[j],muLims[1])
          }
        }
      }
    # } else {
    #   for (i in 1:n) {
    #     familyi <- family[[i]]
    #     if (family == "binomial") {
    #       if (any(mu[i,] < muLims[1] | muLims(2) < mu[i,])) {
    #         for (j in 1:m) {
    #           mu[i,j] <- max(min(mu[i,j],muLims[2]),muLims[1])
    #         }
    #       }
    #     }
    #     if (family == "poisson" || family == "gamma" || 
    #         family == "inverse gaussian") {
    #       if (any(mu[i,] < muLims[1])) {
    #         for (j in 1:m) {
    #           mu[i,j]q() <- max(mu[i,j],muLims[1])
    #         }
    #       }
    #     }
    #   }
    }
    
    # Check stopping conditions
    
    print(max(abs(bvec-bvec.old)))
    if (max(abs(bvec-bvec.old)) < 
        convcrit*max(max(abs(bvec.old))) ) {
      break 
    }
    
  }
  
  #--------------------------------------------------------------------------
  #                    end of GLM iteration loop
  #--------------------------------------------------------------------------
  
  if (iter > iterLim) {
    warning("Iteration limit reached.")
  }
  
  # Sum components of deviance to get the total deviance.
  
  if (is.character(family)) {
    di       <- devFn(mu,y)
    Deviance <- sum((wtvec %*% matrix(1,1,ncurve))*di)
  # } else {
  #   Deviance <- matrix(0,n,ncurve)
  #   for (i in 1:n) {
  #     devFni <- devFn[[i]]
  #     di     <- devFni(mu[i,],y[,i])
  #     Deviance[i,] <- sum((wtvec[i]*matrix(1,1,ncurve))*di)
  #   }
  }
  
  return(list(bvec=bvec, Deviance=Deviance))
  
}

constrain <- function(eta, loBnd, upBnd) {
  eta[eta<loBnd] <- loBnd
  eta[eta>upBnd] <- upBnd
  return(eta)
}

