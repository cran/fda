stepchk <- function(oldstep, cvec, deltac, limwrd, ind,
                    climit=50*c(-rep(1,ncvec), rep(1,ncvec)),
                    active=1:ncvec, dbgwrd) {
  #  stepcheck checks the step size to keep parameters within boundaries
  
  # Last changed 2018 by Jim Ramsay
  
  # define vectors containing lower and upper limits
  
  ncvec   <- length(deltac)
  bot     <- climit[1:ncvec]
  top     <- climit[ncvec+(1:ncvec)]
  
  newstep <- oldstep
  
  #  ensure that step does not go beyond lower limit on parameters
  #  limwrd[2] flags that the lower limit has been hit once
  
  stepi   <- oldstep*deltac
  stepmin <- min(stepi)
  index   <- stepi[active] == stepmin
  if (any(stepi[index] < bot[index]-cvec[index]) &
      any(deltac[index] != 0) )  {
    stepnew <- min((bot[index]-cvec[index])/deltac[index])
    if (dbgwrd) {
      print("Lower limit reached ... new step:")
      cat(c(stepi, round(c(oldstep, stepnew),4)),"\n")
      cat(round(cvec + stepnew*deltac,4),"\n")
    }
    newstep <- stepnew
    if (limwrd[2]) ind <- 1 else limwrd[2] <- TRUE
  } else {
    limwrd[2] <- FALSE
  }
  
  #  check whether upper limit has been reached twice in a row
  
  #  ensure that step does not go beyond upper limit on parameters
  #  limwrd[1] flags that the upper limit has been hit once
  
  stepi   <- oldstep*deltac
  stepmax <- max(stepi)
  index   <- stepi[active] == stepmax
  if (any(stepi[index] > top[index]-cvec[index]) &
      any(deltac[index] != 0) ) {
    stepnew <- min((top[index]-cvec[index])/deltac[index])
    if (dbgwrd) {
      print("Upper limit reached ... new step:")
      cat(c(stepi, round(c(oldstep, stepnew),4)),"\n")
    }
    newstep <- stepnew
    if (limwrd[1]) ind <- 1 else limwrd[1] <- TRUE
  } else {
    limwrd[1] <- FALSE
  }
  
  return(list(newstep, ind, limwrd))
}
