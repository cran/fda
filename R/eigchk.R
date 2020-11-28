eigchk <- function(Cmat) {
    
  #  Last modified 25 August 2020 by Jim Ramsay
  
  #  Cmat for NA's
  
  if (any(is.na(Cmat))) stop("Cmat has NA values.")
  
  #  check Cmat for Cmatmetry
  
  if (max(abs(Cmat-t(Cmat)))/max(abs(Cmat)) > 1e-10) {
    stop('CMAT is not symmetric.')
  } else {
    Cmat <- (Cmat + t(Cmat))/2
  }
  
  #  check Cmat for singularity
  
  eigval <- eigen(Cmat)$values
  ncoef  <- length(eigval)
  if (eigval[ncoef] < 0) {
    neig <- min(length(eigval),10)
    cat("\nSmallest eigenvalues:\n")
    print(eigval[(ncoef-neig+1):ncoef])
    cat("\nLargest  eigenvalues:\n")
    print(eigval[1:neig])
    stop("Negative eigenvalue of coefficient matrix.")
  }
  if (eigval[ncoef] == 0) stop("Zero eigenvalue of coefficient matrix.")
  logcondition <- log10(eigval[1]) - log10(eigval[ncoef])
  if (logcondition > 12) {
    warning("Near singularity in coefficient matrix.")
    cat(paste("\nLog10 Eigenvalues range from\n",
              log10(eigval[ncoef])," to ",log10(eigval[1]),"\n"))
  }
}