varmx.pca.fd <- function(pcafd, nharm = scoresd[2], nx=501)
{
#
#  Apply varimax to the first NHARM components of a pca.fd object.
#  Evaluates the harmonics at NX equally spaced points.
#
#  Returns:
#  An object of class pcafd

#  Note that pcafd is an oldClass type object

#  Last modified 22 August 2006

  if (!(inherits(pcafd, "pca.fd"))) stop(
		"Argument PCAFD is not a pca.fd object.")

  harmfd   <- pcafd$harmonics
  harmcoef <- harmfd$coefs
  coefd    <- dim(harmcoef)
  ndim     <- length(coefd)

  scoresd  <- dim(pcafd$scores)
  if (nharm > scoresd[2]) nharm <- scoresd[2]

  basisobj <- harmfd$basis
  rangex   <- basisobj$rangeval
  x        <- seq(rangex[1], rangex[2], length = nx)
  delta    <- x[2]-x[1]
  harmmat  <- eval.fd(x, harmfd)
  #  If fdmat is a 3-D array, stack into a matrix
  if (ndim == 3) {
     harmmatd <- dim(harmmat)
     dimnames(harmmat) <- NULL
     harmmat  <- aperm(harmmat, c(1, 3, 2))
     dim(harmmat) <- c(harmmatd[1] * harmmatd[3], harmmatd[2])
  }

  #  compute rotation matrix for varimax rotation of harmmat
  rotm <- varmx(harmmat)

  #  rotate coefficients and scores
  if(ndim == 2)
    harmcoef[,1:nharm] <- harmcoef[,1:nharm] %*% rotm
  else
    for(j in (1:coefd[3]))
		harmcoef[,1:nharm,j] <- harmcoef[,1:nharm,j] %*% rotm
  harmscrs <- pcafd$scores[,1:nharm]		
  harmscrs <- harmscrs %*% rotm

  #  compute proportions of variance
  harmvar <- apply(harmscrs^2,2,sum)
  varsum  <- sum(harmvar)
  propvar <- harmvar/varsum

  #  modify pcafd object

  harmfd$coefs    <- harmcoef
  pcafd$harmonics <- harmfd
  pcafd$varprop   <- propvar
  return(pcafd)
}


