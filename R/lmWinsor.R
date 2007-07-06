lmWinsor <- function(formula, data, lower=NULL, upper=NULL, trim=0,
        quantileType=7, subset, weights=NULL, na.action,
        method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
        singular.ok = TRUE, contrasts = NULL, offset=NULL,
        eps=sqrt(.Machine$double.eps), ...)
{
##
## 1.  Identify inputs and outputs
##
  library(quadprog)
#  
  cl <- match.call()
  if(missing(na.action))
    na.action <- get(options('na.action')$na.action) 
  mdly <- mdlx <- formula
  mdly[[3]] <- NULL
  mdlx[[2]] <- NULL
  xNames <- all.vars(mdlx)
  yName <- all.vars(mdly)
  if(length(yName) != 1)
    stop("'formula' must include a single column of 'data'")
  if(as.character(mdly[[2]]) != yName)
    stop("lmWinsor can not accept a formula with a transformed 'y';",
         "  left hand side = ", mdly[[2]], ";  y = ", yName)
##
## 2.  Check 'lower' and 'upper'
##
#  2.1.  Do lower and upper have names?    
  lowNames <- names(lower)
  if(length(lowNames) != length(lower))
    stop("lower must have names")
  hiNames <- names(upper)
  if(length(hiNames) != length(upper))
    stop("upper must have names") 
#  2.2.  Identify numeric columns of 'data' 
  numVars <- sapply(data, is.numeric)
  numV <- names(numVars[numVars])
#  2.3.  Are numeric variables in lower and upper?  
  numLower <- (numV %in% names(lower))
  numUpper <- (numV %in% names(upper))
  nnV <- length(numV)
#  2.4.  Some numeric variables are not in lower;  add   
  {
    if(!all(numLower)){
      Lower <- rep(NA, nnV)
      names(Lower) <- numV
      gotLo <- (numV %in% lowNames) 
      if(any(gotLo)) {
        loGot <- lower[numV[gotLo]]
        Lower[gotLo] <- loGot
      }
      for(v in numV[!numLower])
        Lower[v] <- quantile(data[[v]], trim, na.rm=TRUE, names=FALSE, 
                             type=quantileType)
    }
    else
      Lower <- lower
  }
  {
    if(!all(numUpper)){
      Upper <- rep(NA, nnV)
      names(Upper) <- numV
      gotHi <- (numV %in% hiNames)
      if(any(gotHi)){
        hiGot <- upper[numV[gotHi]]
        Upper[gotHi] <- hiGot
      }
      for(v in numV[!numUpper])
        Upper[v] <- quantile(data[[v]], 1-trim, na.rm=TRUE, names=FALSE, 
                             type=quantileType)
    }
    else
      Upper <- upper
  }
##
## 3.  clipData = data with xNames clipped to (Lower, Upper)
##
  clipData <- data
  for(x. in xNames){
    x.L <- Lower[x.]
    xl <- pmax(data[[x.]], x.L)
#    xl <- pmax(data[[x.]], x.L*(1+3*.Machine$double.eps))
    x.U <- Upper[x.]
    clipData[[x.]] <- pmin(xl, x.U) 
#    clipData[[x.]] <- pmin(xl, x.U *(1-3*.Machine$double.neg.esp))
  }
##
## 4.  fit <- lm(...)
##
  N <- nrow(data)
  if(missing(subset))subset <- 1:N 
#
  cl0 <- as.list(cl)
  cl0[[1]] <- NULL
  cl0$lower <- NULL
  cl0$upper <- NULL
  cl0$trim <- NULL
  cl0$quantileType <- NULL
  cl0$data <- as.name('clipData')
#
  fit <- do.call('lm', cl0)
##
## 5.  Convert to class 'lmWinsor'
##
  fit$call <- cl
  fit$lower <- Lower
  fit$upper <- Upper
  class(fit) <- c("lmWinsor", class(fit))
##
## 6.  all fit$fitted.values %in% (Lower, Upper)[yName]?
##
  y. <- data[[yName]]
  mod <- mean(abs(y.[abs(y.)<Inf]))
  Eps <- {
    if(mod==0) eps
    else eps*mod
  }
  pred <- fit$fitted.value
#  yL <- (Lower[yName]*(1-3*.Machine$double.neg.eps)) 
  yL. <- Lower[yName]
  yL <- yL.-Eps
  yLin <- yL.+Eps 
  yLow <- (pred<yL)
#  yU <- (Upper[yName]*(1+3*.Machine$double.eps)) 
  yU. <- Upper[yName]
  yU <- yU.+Eps 
  yHi <- (yU<pred)
  yUin <- yU.-Eps
# Need yL < yL. < yLin <= yUin < yU. < yU  
# yL == yU?  
  if(yLin>yUin){
    yLin <- yUin <- (yLin+yUin)/2
    yL. <- (yL+yLin)/2
    yU. <- (yU+yUin)/2
  }
#
  if(!any(yLow | yHi)){
    if(!model) fit$model <- NULL
    fit$message <- 'Initial fit in bounds'
    return(fit)
  }
##
## 7.  Else use quadratic progamming to minimize the
##     Winsorized sum of squares of residuals 
##
  out <- matrix(FALSE, N, 2, dimnames=list(NULL,
                       c('out', 'high')))
  coef0 <- coef(fit)
  k <- sum(!is.na(coef0))
  extraStats <- c("SSEraw", "SSEclipped",
      "nLoOut", "nLo.", "nIn", "nHi.", "nHiOut")   
  coefiter <- matrix(NA, N+1, k+7, dimnames=list(NULL,
        c(names(coef0[!is.na(coef0)]), extraStats)) )
  coefiter[1, 1:k] <- coef0[!is.na(coef0)]
  coefiter[1, "SSEraw"] <- sum((y.-pred)^2)
  predW <- pmax(yL., pmin(yU., pred)) 
  coefiter[1, "SSEclipped"] <- sum((y.-predW)^2)
  coefiter[1, "nLoOut"] <- sum(yLow)
  coefiter[1, "nLo."] <- sum(!yLow & (pred<yLin))
  coefiter[1, "nHiOut"] <- sum(yHi)
  coefiter[1, "nHi."] <- sum(!yHi & (yUin<=pred))
  coefiter[1, "nIn"] <- (N-sum(coefiter[1,
           c("nLoOut", "nLo.", "nHi.", "nHiOut") ])) 
#
  templimits <- matrix(NA, N, 3, dimnames=list(NULL,
                       c("newConstraint", "lower", "upper") ) )
#  
  X <- model.matrix(fit)[, !is.na(coef0), drop=FALSE]
  yMin <- min(y.)
  yMax <- max(y.) 
#   
  for(i in 1:N){
#   The standard exit from this loop is via 'break'
#   afer all constraints are satisfied.
#   7.1.  Find prediction farthest out 
    dLow <- (yL-pred)
    dHi <- (pred-yU)
    dLow. <- dLow[!out[, 1]]
    dHi. <- dHi[!out[, 1]] 
    {
      if(max(dLow.) > max(dHi.)){
        lows <- (which(!out[, 1])[max(dLow.)==dLow.])
        lowest <- lows[y.[lows]==min(y.[lows])]
        out[lowest, ] <- c(TRUE, FALSE)
        templimits[i, "newConstraint"] <- lowest
        yMin <- min(yLin, pred[pred>pred[lowest]])
      }
      else{
        highs <- (which(!out[, 1])[max(dHi.)==dHi.])
        highest <- highs[y.[highs]==max(y.[highs])]      
        out[highest, ] <- TRUE
        templimits[i, "newConstraint"] <- highest
        yMax <- max(yUin, pred[pred<pred[highest]]) 
      }
    }
#   7.2.  Use QP ...
#    Dmat. <- crossprod(X[!out[, 1], , drop=FALSE])
    qrX <- qr(X[!out[, 1], , drop=FALSE])
    qR <- qr.R(qrX)
    signR <- sign(diag(qR))
    qR. <- (signR * qR)
#    qrXdiag <- diag(qR.) 
    qrXrng <- range(diag(qR.))
    sing <- ((dim(qR)[1] < k) || (qrXrng[1] < eps * qrXrng[2]))
    if(sing){
      fit$message <- 'Iteration terminated by a singular quadratic program'
      fit$qr <- {
        if(i<2) qrX else qrX.old
      }
      break
    }
    Dmat. <- solve(qR.) 
#
    dvec. <- as.numeric(crossprod(X[!out[, 1], , drop=FALSE],
                                  y.[!out[, 1]])) 
#
    outL <- (out[, 1] & !out[, 2])
    outH <- (out[, 1] & out[, 2]) 
    AmatL <- X[outL, , drop=FALSE]
    AmatH <- X[outH, , drop=FALSE]
    Amat. <- rbind(-AmatL, AmatH)
#
#    maxPredLo <- max(pred[outL], -Inf)
#    minPredHi <- min(pred[outH], Inf)
#    yMin <- min(y.[y.>maxPredLo], yL)+Eps
#
#    yMax <- max(y.[y.<minPredHi], yU)-Eps
#
    bvec. <- c(rep(-yMin, nrow(AmatL)), rep(yMax, nrow(AmatH)))
#
    templimits[i, c("lower", "upper")] <- c(yMin, yMax)
#
#    QPi <- solve.QP(Dmat=Dmat., dvec=dvec., Amat=Amat., bvec=bvec.)
    env <- new.env()
    assign("D", Dmat., envir=env)
    assign("d", dvec., envir=env)
    assign("A", t(Amat.), envir=env)
    assign("b", bvec., envir=env)
    QPlist <- list(Dmat=quote(D), dvec=quote(d),
                   Amat=quote(A), bvec=quote(b),
                   factorized=TRUE)
    QPi <- do.call('solve.QP', QPlist, envir=env)
#    QPi <- solve.QP(Dmat., as.numeric(dvec), Amat, bvec)
#   7.3.  Are unconstrained predictions inside?
#    pred.old <- pred
    i1 <- i+1
    coefiter[i1, 1:k] <- QPi$solution 
    pred <- X %*% QPi$solution 
#    pred <- pmax(yL, pmin(yU, pred0))
    yLow <- (pred<yL)
    yHi <- (yU<pred)
    coefiter[i1, "SSEraw"] <- sum((y.-pred)^2)
    predW <- pmax(yL., pmin(yU., pred)) 
    coefiter[i1, "SSEclipped"] <- sum((y.-predW)^2)
    coefiter[i1, "nLoOut"] <- sum(yLow)
    coefiter[i1, "nLo."] <- sum(!yLow & (pred<yLin))
    coefiter[i1, "nHiOut"] <- sum(yHi)
    coefiter[i1, "nHi."] <- sum(!yHi & (yUin<=pred))
    coefiter[i1, "nIn"] <- (N-sum(coefiter[i1,
           c("nLoOut", "nLo.", "nHi.", "nHiOut") ])) 
#    
    out2 <- ((pred < yL) | (yU < pred))
    if(!any(out2[!out[, 1]])) {
      fit$message <- 'QP iterations successful'
      fit$qr <- qrX 
      break
    }
    qrX.old <- qrX 
#   7.4.  Find the next extreme prediction ... -> 7.1 
  }
##
## 8.  Modify fit as appropriate
##
  coef0[!is.na(coef0)] <- QPi$solution
  fit$coefficients <- coef0
#
  predy <- pmax(yL, pmin(yU, pred))
  fit$residuals <- (y.-predy)
  fit$fitted.values <- predy
  fit$coefIter <- coefiter[!is.na(coefiter[, 1]), , drop=FALSE]
  fit$tempLimits <- templimits[!is.na(templimits[, 1]), , drop=FALSE]
#
  fit
}
