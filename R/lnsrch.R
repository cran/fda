lnsrch <- function (xold, fold, g, p, func, dataList, stpmax, 
                    itermax=20, TOLX=1e-10, dbglev=0) {
    n     <- length(xold)
    check <- FALSE
    f2    <- 0
    alam2 <- 0
    ALF   <- 1e-4
    psum  <- sqrt(sum(p^2))
    #  scale if attempted step is too big
    if (psum > stpmax) {
        p <- p*(stpmax/psum)
    }
    slope <- sum(g*p)
    if (dbglev > 1) {
      cat("\n")
      cat("      ")
      cat(0)
      cat("      ")
      cat(round(slope,5))
      cat("      ")
      cat(round(fold,5))
    } 
    if (slope >= 0) {
        stop('Initial slope not negative.')
    }
    # compute lambdamin
    test <- 0
    for (i in 1:n) {
        temp <- abs(p[i])/max(abs(xold[i]),1)
        if (temp > test) {
            test <- temp
        }
    }
    alamin <- TOLX/test
    #  always try full Newton step first
    alam   <- 1
    #  start of iteration loop
    iter <- 0
    while (iter <= itermax) {
        iter <- iter + 1
        x <- xold + alam*p
        #  -------------  function evaluation  -----------
        result <- func(x, dataList)
        f <- result$PENSSE
        g <- result$DPENSSE
        if (dbglev > 1) {
          cat("\n")
          cat("      ")
          cat(iter)
          cat("      ")
          cat(round(slope,5))
              cat("      ")
              cat(round(fold,5))
        }
        #  -----------------------------------------------
        #  convergence on x change.
        if (alam < alamin) {
            x     <- xold
            check <- TRUE
            return(list(x=x, check=check))
        } else {
            #  sufficient function decrease
            if (f <= fold + ALF*alam*slope) {
                return(list(x=x, check=check))
            }
            #  backtrack
            if (alam == 1) {
                #  first time
                tmplam <- -slope/(2*(f-fold-slope))
            } else {
                #  subsequent backtracks
                rhs1 <- f  - fold - alam *slope
                rhs2 <- f2 - fold - alam2*slope
                a <- (rhs1/alam^2 - rhs2/alam2^2)/(alam-alam2)
                b <- (-alam2*rhs1/alam^2 + alam*rhs2/(alam*alam2))/(alam-alam2)
                if (a == 0) {
                    tmplam <- -slope/(2*b)
                } else {
                    disc <- b^2 - 3*a*slope
                    if (disc < 0) {
                        tmplam <- 0.5*alam
                    } else {
                        if (b <= 0) {
                            tmplam <- (-b+sqrt(disc))/(3*a)
                        } else {
                            tmplam <- -slope/(b+sqrt(disc))
                        }
                    }
                    if (tmplam > 0.5*alam) {
                        tmplam <- 0.5*alam
                    }
                }
            }
            alam2 <- alam
            f2    <- f
            #  lambda > 0.1 lambda1
            alam <- max(tmplam, 0.1*alam)
        }
        #  try again
    }
    
    return(list(x=x, check=check))
    
}

