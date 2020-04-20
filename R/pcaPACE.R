pcaPACE <- function(covestimate, nharm = 2, harmfdPar=NULL, cross = TRUE)
{
  #  Carry out a functional PCA with regularization from the estimate of the covariance surface
  #  Arguments:
  #  COVESTIMATE  ... list of length 2 with named entries cov.estimate and meanfd
  #  NHARM     ... Number of principal components or harmonics to be kept
  #  HARMFDPAR ... Functional parameter object for the harmonics
  #  CROSS     ... If TRUE, ........................
  #
  #  Returns:  An object PCAFD of class "pca.PACE.fd" with these named entries:
  #  harmonics  ... A functional data object for the harmonics or eigenfunctions
  #  values     ... The complete set of eigenvalues
  #  scores     ... TO DO ===================================
  #  varprop    ... A vector giving the proportion of variance explained
  #                 by each eigenfunction
  #
  
  #  set up HARMBASIS
  
  harmbasis <- harmfdPar$fd$basis
  nhbasis   <- harmbasis$nbasis
  
  #  set up LFDOBJ and LAMBDA
  
  Lfdobj <- harmfdPar$Lfd
  lambda <- harmfdPar$lambda
  
  
  #  set up cross product Lmat for harmonic basis,
  #  roughness penalty matrix Rmat, and
  #  penalized cross product matrix Lmat
  
  Lmat <- eval.penalty(harmbasis, 0)
  if (lambda > 0) {
    Rmat <- eval.penalty(harmbasis, Lfdobj)
    Lmat <- Lmat + lambda * Rmat
  }
  Lmat <- (Lmat + t(Lmat))/2
  
  #  compute the Choleski factor Mmat of Lmat
  
  Mmat    <- chol(Lmat)
  Mmatinv <- solve(Mmat)
  
  #  set up cross product and penalty matrices
  
  if(class(covestimate$cov.estimate) == "bifd"){
    Wmat = covestimate$cov.estimate$coefs
    nvar = 1
    basisobj = covestimate$cov.estimate$sbasis
  }else{
    k=0
    l = length(covestimate$cov.estimate)
    while(l>0){
      k = k + 1
      l = l - k
    }
    nvar= k
    basisobj = covestimate$cov.estimate[[1]]$sbasis
    nbasis = basisobj$nbasis
    
    if(!cross){
      diag = cumsum(c(1,k - 0:(k-2)))
      Wmat = lapply(covestimate$cov.estimate[diag], function(x) x$coefs)
    }else{
      Wmat = matrix(0,nvar*nbasis,nvar*nbasis)
      r=0
      for(i in 1:nvar){
        indexi <- 1:nbasis + (i - 1) * nbasis
        for(j in 1:nvar){
          indexj <- 1:nbasis + (j - 1) * nbasis
          if(j>=i){
            r=r+1
            Wmat[indexi,indexj] = covestimate$cov.estimate[[r]]$coefs
            Wmat[indexj,indexi] = covestimate$cov.estimate[[r]]$coefs
          }
        }
      }
    }
    
    
    
  }
  
  Jmat = inprod(harmbasis, basisobj)
  MIJW = crossprod(Mmatinv,Jmat)
  
  #  set up matrix for eigenanalysis
  nbasis = basisobj$nbasis
  if(nvar == 1) {
    Cmat = MIJW %*% Wmat %*% t(MIJW)
  } else {
    if(!cross){
      Cmat = lapply(Wmat, function(x) MIJW %*% x %*% t(MIJW) )
      Cmat = lapply(Cmat, function(x) (x + t(x))/2)
    }else{
      Cmat = matrix(0,nvar*nhbasis,nvar*nhbasis)
      for(i in 1:nvar){
        indexi <- 1:nbasis + (i - 1) * nbasis
        for(j in 1:nvar){
          indexj <- 1:nbasis + (j - 1) * nbasis
          Cmat[indexi,indexj] = MIJW %*% Wmat[indexi,indexj] %*% t(MIJW)
          Cmat[indexj,indexi] = MIJW %*% Wmat[indexj,indexi] %*% t(MIJW)
        }
      }
    }
    
    
  }

  
  #  eigenalysis
  
  if (nvar == 1 | cross){
    Cmat    <- (Cmat + t(Cmat))/2
    result  <- eigen(Cmat)
    eigvalc <- result$values[1:nharm]
    if(nvar>1) eigvalc <- t(replicate(nvar, eigvalc))
    eigvecc <- as.matrix(result$vectors[, 1:nharm])
    sumvecc <- apply(eigvecc, 2, sum)
    eigvecc[,sumvecc < 0] <-  - eigvecc[, sumvecc < 0]
    varprop <- eigvalc[1:nharm]/sum(eigvalc) 
  } else {
    Cmat  = lapply(Cmat, function(x) (x + t(x))/2) 
    result = lapply(Cmat,eigen)
    eigval = lapply(result, function(x) x$values[1:nharm])
    eigvalc = do.call(rbind, eigval)
    eigvecc = lapply(result, function(x) as.matrix(x$vectors[, 1:nharm]))
    
    sumvecc = lapply(eigvecc, function(x) apply(x,2,sum))
    for(r in 1:length(eigvecc)){
      eigvecc[[r]][,sumvecc[[r]]<0] <- - eigvecc[[r]][,sumvecc[[r]]<0]
    }
    varprop <- eigvalc[1:nharm]/sum(eigvalc)
  }
  
  #  set up harmfd
  
  if (nvar == 1) {
    harmcoef <- Mmatinv %*% eigvecc
  } else {
    harmcoef <- array(0, c(nbasis, nharm, nvar))
    if(!cross){
      for (j in 1:nvar) {
        harmcoef[,  , j] <- Mmatinv %*% eigvecc[[j]]
      }
    }else{
      for (j in 1:nvar) {
        index <- 1:nbasis + (j - 1) * nbasis
        harmcoef[,  , j] <- Mmatinv %*% eigvecc[index,  ]
      }
    }
    
    
  }
  harmnames <- rep("", nharm)
  for(i in 1:nharm)
    harmnames[i] <- paste("PC", i, sep = "")
  harmfd   <- fd(harmcoef, harmbasis)
  
  #  set up the object pcafd of the pca.fd class containing the results
  
  scores = NULL
  pcafd        <- list(harmfd, eigvalc, scores,varprop, covestimate$meanfd)
  class(pcafd) <- "pca.fd"
  names(pcafd) <- c("harmonics", "values","scores" ,"varprop", "meanfd")
  
  return(pcafd)
}

