covPACE <- function(data,rng , time, meanfd, basis, lambda, Lfdobj){
  # does a bivariate smoothing for estimating the covariance surface for data that has not yet been smoothed
  #
  #   Arguments:
  #   DATA .....  a matrix object or list -- If the set is supplied as a matrix object, 
  #               the rows must correspond to argument values and columns to replications, 
  #               and it will be assumed that there is only one variable per observation.  
  #               If y is a three-dimensional array, the first dimension corresponds to  
  #               argument values, the second to replications, and the third to variables 
  #               within replications. -- If it is a list, each element must be a matrix
  #               object, the rows correspond to argument values per individual. First 
  #               column corresponds to time points and the following columns to argument  
  #               values per variable.
  #   RNG .....   a vector of length 2 defining a restricted range where the data was observed
  #   TIME .....  Array with time points where data was taken. length(time) == dim(data)[1]
  #   MEANFD .... Fd object corresponding to the mean function of the data
  #   BASIS ..... basisfd object for smoothing the covariate function
  #   LAMBDA .... a nonnegative real number specifying the amount of smoothing to be applied to the estimated
  #               functional parameter
  #   Lfdobj .... linear differential operator object for smoothing penalty of the estimated functional parameter
  #
  # Returns a list with the two named entries cov.estimate and meanfd
  
  if(is.list(data)){
    datalist = data 
    s = length(data)
  }else{
    datalist = sparse.list(data,time)
    s = dim(data)[2]
  }
  
  indexes = lapply(datalist, function(x) which(time %in% x[,1])) #sampling points for each subject
  nvar = dim(datalist[[1]])[2] - 1
  
  #centering the data
  mean.point = matrix(eval.fd(time, meanfd),nrow = nrow(eval.fd(time, meanfd)),ncol=nvar)
  data.center = lapply(datalist, function(x){
    y = x[,-1] - mean.point[which(time %in% x[,1]),] 
    return(y)
  })
  data.center = lapply(data.center, as.matrix)

  #create kronecher product matrix for basis evaluation
  phi.eval = eval.basis(time,basis)
  phi.kro = kronecker(phi.eval,phi.eval)
  m = length(time)
  phi.kro.sub = lapply(m*(1:m)-m+1, function(x) return(phi.kro[x:(x+m-1),])) 
  phi.kro.sub = lapply(phi.kro.sub,as.matrix)
  
  #point covariance estimate 
  
  
  if(nvar > 1){
    out = NULL
    for(i in 1:(nvar-1)){
      out = c(out,i:(nvar-1)*nvar+i)
    }
    
    pairs = (expand.grid(1:nvar,1:nvar))[-out,]
    z = list()  
    for(q in 1:nrow(pairs)){
      y = pairs[q,]
      z[[q]] = lapply(data.center, function(x) as.vector(x[,pairs[q,1]] %*% t(x[,pairs[q,2]]))[-(1:nrow(x) + nrow(x) *(0:(nrow(x)-1)))])
    }
  }else{
  z = lapply(data.center, function(x) as.vector(x %*% t(x))[-(1:nrow(x) + nrow(x) *(0:(nrow(x)-1)))])
  }
  
  
  #Calculate t(phi)%*%phi and t(phi)%*%y
  
  phi.norm = 0 # t(phi)%*%phi -- nbasis^2 x nbasis^2 matrix 
  if(nvar == 1) phi.y = 0 #t(phi)%*%y -- nbasis^2 x 1
  if(nvar > 1) phiy.list = as.list(replicate(length(z),0))
  
  for (j in 1:s){
    l = length(indexes[[j]])
    phij = do.call(rbind,lapply(phi.kro.sub[indexes[[j]]], function(x) x[indexes[[j]],]))
    phij = phij[-(1:l + l *(0:(l-1))),] #take out the diagonal 
    phi.norm = phi.norm + t(phij)%*%phij
    
    if(nvar > 1){
      phiy.list = lapply(1:length(z), function(x) phiy.list[[x]] + t(phij)%*%z[[x]][[j]])
    }else{
      phi.y = phi.y + t(phij)%*%z[[j]]
    }
  }
  
  if(nvar > 1) phi.y = do.call(cbind,phiy.list)
  
  
  #Penalty matrix
  Pl = eval.penalty(basis,Lfdobj,rng)
  Po = inprod(basis,basis,rng = rng)
  P = kronecker(t(Pl),Po) + kronecker(Po,Pl)
  
  
  c.estimate = solve(phi.norm+lambda*P) %*% (phi.y)
  
  if(nvar > 1){
    c.list = lapply(1:length(z), function(x) matrix(c.estimate[,x],nrow = sqrt(nrow(c.estimate)), ncol = sqrt(nrow(c.estimate))) )
    cov.estimate = lapply(c.list, function(x) bifd(x,basis,basis) )
  }else{
    c = matrix(c.estimate,nrow = sqrt(length(c.estimate)), ncol = sqrt(length(c.estimate)))
    cov.estimate = bifd(c,basis,basis)
  }
  
  covPACE = list(cov.estimate,meanfd)
  names(covPACE) = c("cov.estimate","meanfd")
  return(covPACE)
}



