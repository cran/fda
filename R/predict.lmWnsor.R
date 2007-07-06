predict.lmWinsor <- function(object, newdata, se.fit = FALSE,
     scale = NULL, df = Inf, 
     interval = c("none", "confidence", "prediction"),
     level = 0.95, type = c("response", "terms"), 
     terms = NULL, na.action = na.pass,
     pred.var = res.var/weights, weights = 1, ...){
##
## 0.  Set up
##  
  cl <- match.call()
# Define res.var so R CMD check won't complain
  res.var <- 0
##
## 1.  Identify inputs and outputs
##
  mdly <- mdlx <- formula(object) 
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
## 2.  Check limits 
##
  Lower <- object$lower
  got.xL <- (xNames %in% names(Lower))
  if(any(!got.xL))
    stop("No lower limit for ",
         paste(names(got.xL[!got.xL]), collapse=", "))
  got.yL <- (yName %in% names(Lower))
  if(!got.yL)
    stop("No lower limit for ", yName)
#  
  Upper <- object$upper
  got.xU <- (xNames %in% names(Upper))
  if(any(!got.xU))
    stop("No upper limit for ",
         paste(names(got.xU[!got.xU]), collapse=", "))
  got.yU <- (yName %in% names(Upper))
  if(!got.yU)
    stop("No upper limit for ", yName)
##
## 3.  Clip newdata[, xNames] and predict 
##  
  {
    if(missing(newdata))
      pred <- object$fitted.values
    else{
      newDat <- newdata
      for(x. in xNames){
        xl <- pmax(newdata[[x.]], Lower[x.])
        newDat[[x.]] <- pmin(xl, Upper[x.])
      }
      cl0 <- as.list(cl)
      cl0[[1]] <- NULL
      cl0$newdata <- as.name('newDat')       
      pred <- do.call('predict.lm', cl0)      
    }
  }
##
## 4.  Clip pred
##
  yl <- pmax(pred, Lower[yName])
  y. <- pmin(yl, Upper[yName])
##
## 5.  done
##
  y.
}
