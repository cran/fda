yListCheck = function(yList, nvar) {
  
  #  Last modified 3 August 2017
  
  if (!is.list(yList)) {
    stop("YLIST is not a list.")
  }
  errwrd = FALSE
  
  nvec    = rep(0,nvar)
  dataWrd = rep(FALSE,nvar)
  ydim    = rep(0,nvar)
  nobs    = 0
  #  find number of replications for (first non-empty cell
  for (ivar in 1:nvar) {
    if (is.list(yList[[ivar]])) {
      yListi = yList[[ivar]]
      nrep = dim(as.matrix(yListi$y))[2]
      break
    }
  }
  #  loop through variables
  for (ivar in 1:nvar) {
    if (is.list(yList[[ivar]]) && !is.null(yList[[ivar]]$y)) {
      dataWrd[ivar] = TRUE
      yListi = yList[[ivar]]
      if (is.null(yListi$argvals)) {
        warning(paste("ARGVALS is not a member for (YLIST[[", ivar,"]]."))
        errwrd = TRUE
      }
      ni = length(yListi$argvals)
      nvec[ivar] = ni
      ydimi = dim(as.matrix(yListi$y))
      if (length(ydimi) > 2) {
        warning(paste("More than two dimensions for (y in YLIST[[",
                      ivar,"]]."))
        errwrd = TRUE
      } else {
        ydim[ivar] = ydimi[1]
      }
      #  set up and check NREP
      nrepi = ydimi[2]
      if (nrepi != nrep) {
        warning("Second dimensions of YList.y are not equal.")
        errwrd = TRUE
      }
      nobs = nobs + 1
      if (ni != ydimi[1]) {
        print(c(ni,ydimi[1]))
        warning(paste("Length of ARGVALS and first dimension of Y",
                      "are not equal."))
        errwrd = TRUE
      }
    } else {
      dataWrd[ivar] = FALSE
    }
  }
  
  if (nobs == 0) {
    warning("No variables have observations.")
    errwrd = TRUE
  }
  
  if (errwrd) {
    stop("One or more terminal stop encountered in YLIST.")
  }
  
  return(list(nrep = nrep, nvec = nvec, dataWrd = dataWrd))
  
}
