funcint <- function(func,cvec, basisobj, nderiv=0, JMAX=15, EPS=1e-7)
{
  
  #  computes the definite integral of a function defined in terms of 
  #  a functional data object with basis BASISOBJ and coefficient vector CVEC.
  
  #  FUNC     ... the name of a function
  #  CVEC     ..  a vector of coefficients defining the functional data
  #               object required to compute the value of the function.
  #  BASISOBJ ... a functional data basis object
  #  NDERIV   ... a non-negative integer defining a derivative to be used.
  #  JMAX   maximum number of allowable iterations
  #  EPS    convergence criterion for relative stop
  
  #  Return: the value of the function
  
  #  Last modified 15 June 2020 by Jim Ramsay
  
  #  set iter
  
  iter <- 0
  
  # The default case, no multiplicities.
  
  rng    <- basisobj$rangeval
  
  nbasis <- basisobj$nbasis
  
  nrep   <- dim(basisobj$coef)[2]
  
  inprodvec <- matrix(0,nrep,1)
  
  #  set up first iteration
  
  iter  <- 1
  width <- rng[2] - rng[1]
  JMAXP <- JMAX + 1
  h     <- rep(1,JMAXP)
  h[2]  <- 0.25
  s <- matrix(0,JMAXP,nrep)
  sdim <- length(dim(s))
  #  the first iteration uses just the endpoints
  fdobj <- fd(cvec, basisobj)
  fx    <- func(eval.fd(rng, fdobj, nderiv))
  #  multiply by values of weight function if necessary
  s[1,] <- width*apply(fx,2,sum)/2
  tnm  <- 0.5
  
  #  now iterate to convergence
  
  for (iter in 2:JMAX) {
    tnm <- tnm*2
    if (iter == 2) {
      x <- mean(rng)
    } else {
      del <- width/tnm
      x   <- seq(rng[1]+del/2, rng[2]-del/2, del)
    }
    fx <- func(eval.fd(x, fdobj, nderiv))
    chs <- width*apply(fx,2,sum)/tnm
    s[iter,] <- (s[iter-1,] + chs)/2
    if (iter >= 5) {
      ind <- (iter-4):iter
      ya <- s[ind,,]
      ya <- matrix(ya,5,nrep)
      xa <- h[ind]
      absxa <- abs(xa)
      absxamin <- min(absxa)
      ns <- min((1:length(absxa))[absxa == absxamin])
      cs <- ya
      ds <- ya
      y  <- ya[ns,,]
      ns <- ns - 1
      for (m in 1:4) {
        for (i in 1:(5-m)) {
          ho      <- xa[i]
          hp      <- xa[i+m]
          w       <- (cs[i+1,] - ds[i,])/(ho - hp)
          ds[i,] <- hp*w
          cs[i,] <- ho*w
        }
        if (2*ns < 5-m) {
          dy <- cs[ns+1,]
        } else {
          dy <- ds[ns,]
          ns <- ns - 1
        }
        y <- y + dy
      }
      ss     <- y
      errval <- max(abs(dy))
      ssqval <- max(abs(ss))
      if (all(ssqval > 0)) {
        crit <- errval/ssqval
      } else {
        crit <- errval
      }
      if (crit < EPS && iter >= 5) break
    }
    s[iter+1,] <- s[iter,]
    h[iter+1]   <- 0.25*h[iter]
    if (iter == JMAX) warning("Failure to converge.")
  }
  inprodvec <- inprodvec + ss
  
}


