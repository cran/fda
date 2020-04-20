reconsCurves <- function(data, PC){
  #Reconstruct data curves using functional principal components
  #
  # Arguments:
  #
  # DATA ...... an set of values of curves at discrete sampling points or 
  #             argument values. If the set is supplied as a matrix object, 
  #             the rows must correspond to argument values and columns to 
  #             replications, and it will be assumed that there is only one 
  #             variable per observation. If data is a three-dimensional array, 
  #             the first dimension corresponds to argument values, the second 
  #             to replications, and the third to variables within replications.
  # PC ......   an object of class "pca.fd"
  #
  # Returns a functional data object (i.e., having class "fd")
  ndim = length(dim(data))
  
  if(ndim == 3){
    nvar = dim(data)[3]
    nbasish = PC$harmonics$basis$nbasis
    datarecon = 0
    data.coefs = array(NA, dim = c(nbasish,ncol(data),nvar))
    for(i in 1:ncol(PC$harmonics$coefs)){
      for(l in 1:nvar){
        data.coefs[,,l] = t(replicate(nbasish,PC$scores[,i,l]))*PC$harmonics$coefs[,i,l]
      }
      datarecon = datarecon + fd(data.coefs,PC$harmonics$basis)
    }
    
    
    reconlist = list()
    for(p in 1:nvar){
      meanfdaux = PC$meanfd
      meanfdaux$coefs = replicate(ncol(data),PC$meanfd$coefs[,,p])
      reconaux = datarecon
      reconaux$coefs = as.matrix(datarecon$coefs[,,p])
      reconlist[[p]] = meanfdaux + reconaux
    }
    
    finalcoef = array(NA, dim = c(nrow(reconlist[[1]]$coefs), ncol(data),nvar))
    for(n in 1:nvar){
      finalcoef[,,n] = reconlist[[n]]$coefs
    }
    
    data.recons = reconlist[[1]]
    data.recons$coefs = finalcoef
    data.recons$fdnames$values = data.recons$fdnames$reps
    
  }else{
    datarecon = 0
    for(i in 1:ncol(PC$harmonics$coefs)){
      coefs = t(replicate(PC$harmonics$basis$nbasis,PC$scores[,i]))*PC$harmonics$coefs[,i]
      if(dim(data)[2]==1) 
      coefs = replicate(PC$harmonics$basis$nbasis,PC$scores[,i])*PC$harmonics$coefs[,i]
      datarecon = datarecon + fd(coefs,PC$harmonics$basis)
    }
    PC$meanfd$coefs = replicate(ncol(datarecon$coefs),as.vector(PC$meanfd$coefs))
    data.recons = PC$meanfd+datarecon
    data.recons$fdnames$values = "values"
  }
  
  
  return(data.recons)
}
