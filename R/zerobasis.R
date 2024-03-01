zerobasis <- function(k) {
# ZEROBASIS constructes a K by K-1 matrix that maps an unrestricted matrix B with K - 1 rows by 
#  the linear transformation ZEROBASIS %*% B = C into the subspace of matrices with K rows having #  column sums equal to zero.  
#  The matrix has orthonormal columns, so that crossprod(ZEROBASIS) is the identity matrix
#  of order K - 1.

  tk <- 0:(k-1) + 0.5
  fbasis     <- create.fourier.basis(k,k)
  fbasmat    <- eval.basis(tk, fbasis)
  if (k > 2) fbasmat <- fbasmat[,2:k] else fbasmat <- matrix(c(1,-1),2,1)
  fbasnorm   <- sqrt(apply(fbasmat^2,2,sum))
  zerobasmat <- fbasmat/outer(rep(1,k),fbasnorm)
  return(zerobasmat)
}