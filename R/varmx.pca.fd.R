varmx.pca.fd <- function(pcafd, nharm = scoresd[2], nx=101)
{
#
#  Apply varimax to the first NHARM components of a pca.fd object.
#  Evaluates the harmonics at NX equally spaced points.
#
#  Returns:
#  An object of class pcafd

#  Last modified 23 January 2003

  if (!(inherits(pcafd, "pca.fd"))) stop("Argument PCAFD is not a pca.fd object.")

  harmfd   <- pcafd[[1]]
  harmcoef <- getcoef(harmfd)
  coefd    <- dim(harmcoef)
  ndim     <- length(coefd)

  scoresd  <- dim(pcafd$scores)
  if (nharm > scoresd[2]) nharm <- scoresd[2]

  basisfd  <- getbasis(harmfd)
  rangex   <- basisfd$rangeval
  x        <- seq(rangex[1], rangex[2], length = nx)
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
  if(ndim == 2) {
    harmcoef <- harmcoef[,1:nharm]
    harmcoef <- harmcoef %*% rotm
  } else {
    harmcoef <- harmcoef[,1:nharm,]
    for(j in (1:coefd[3])) harmcoef[,,j] <- harmcoef[,,j] %*% rotm
  }
  pcafd$scores  <- pcafd$scores[,1:nharm] %*% rotm
  #  modify pcafd object
  harmfd[[1]] <- harmcoef
  pcafd[[1]]  <- harmfd
  pcafd$varprop <- apply(pcafd$scores^2, 2, mean)/sum(pcafd$values)
  setOldClass("pca.fd")  
  oldClass(pcafd) <- "pca.fd"
  return(pcafd)
}

varmx <- function(amat) {

  #  Does a VARIMAX rotation of a principal components solution

  #  Arguments:
  #  AMAT  ...  N by K matrix of harmonic values

  #  Returns:
  #  ROTM  ...  Rotation matrixed loadings

  n <- nrow(amat)
  k <- ncol(amat)
  rotm <- diag(k)
  if (k == 1) return(rotm)

  eps  <- 0.0011
  ccns <- 0.7071068

  varold <- 0
  varnow <- sum(apply(amat^2, 2, var))

  iter <- 0
  while (abs(varnow - varold) > 1e-7 && iter <= 50) {
    iter  <- iter + 1
    for (j in 1:(k-1)) for (l in (j+1):k) {
      avecj  <- amat[,j]
      avecl  <- amat[,l]
      uvec   <- avecj^2 - avecl^2
      tvec   <- 2*avecj*avecl
      aa <- sum(uvec)
      bb <- sum(tvec)
      cc <- sum(uvec^2 - tvec^2)
      dd <- 2*sum(uvec*tvec)
      tval <- dd - 2*aa*bb/n
      bval <- cc - (aa^2 - bb^2)/n

      if (tval == bval) {
        sin4t <- ccns
        cos4t <- ccns
      }

      if (tval < bval) {
        tan4t <- abs(tval/bval)
        if (tan4t >= eps) {
          cos4t <- 1/sqrt(1+tan4t^2)
          sin4t <- tan4t*cos4t
        } else {
          if (bval < 0) {
            sin4t <- ccns
            cos4t <- ccns
          } else {
            sin4t <- 0
            cos4t <- 1
          }
        }
      }

      if (tval > bval) {
        ctn4t <- abs(tval/bval)
        if (ctn4t >= eps) {
          sin4t <- 1/sqrt(1+ctn4t^2)
          cos4t <- ctn4t*sin4t
        } else {
          sin4t <- 1
          cos4t <- 0
        }
      }

      cos2t <- sqrt((1+cos4t)/2)
      sin2t <- sin4t/(2*cos2t)
      cost  <- sqrt((1+cos2t)/2)
      sint  <- sin2t/(2*cost)
      if (bval > 0) {
        cosp <- cost
        sinp <- sint
      } else {
        cosp <- ccns*(cost + sint)
        sinp <- ccns*abs(cost - sint)
      }
      if (tval <= 0) sinp <- -sinp

      amat[,j] <-  avecj*cosp + avecl*sinp
      amat[,l] <- -avecj*sinp + avecl*cosp
      rvecj    <- rotm[,j]
      rvecl    <- rotm[,l]
      rotm[,j] <-  rvecj * cosp + rvecl * sinp
      rotm[,l] <- -rvecj * sinp + rvecl * cosp

    }
    varold <- varnow
    varnow <- sum(apply(amat^2,2,var))
  }

  return( rotm )
}
