smooth.sparse.mean <- function(data, time ,rng=c(0, 1), type = "" , nbasis = NULL, knots = NULL, norder = NULL, lambda = NULL){
  
  #   Arguments:
  #   DATA        a matrix object or list -- If the set is supplied as a matrix object, 
  #               the rows must correspond to argument values and columns to replications, 
  #               and it will be assumed that there is only one variable per observation.  
  #               If y is a three-dimensional array, the first dimension corresponds to  
  #               argument values, the second to replications, and the third to variables 
  #               within replications. -- If it is a list, each element must be a matrix
  #               object, the rows correspond to argument values per individual. First 
  #               column corresponds to time points and followin columns to argument values 
  #               per variable.
  #   TIME        Array with time points where data was taken. length(time) == ncol(data)
  #   RNG         an array of length 2 containing the lower and upper
  #               boundaries for the rangeval of argument values
  #   TYPE        Type of basisfd for smoothing the mean estimate function
  #   NBASIS      An integer variable specifying the number of basis functions
  #   KNOTS       a vector specifying the break points if type = "bspline"
  #   NORDER      an integer specifying the order of b-splines if type = "bspline"
  #   LAMBDA      a nonnegative real number specifying the amount of smoothing to be applied 
  #               to the estimated functional parameter
  
  if(type == "bspline"){
    nbasis = length(knots) + norder - 2
    basis = create.bspline.basis(rng,nbasis,norder)
  }else if (type == "fourier"){
    basis = create.fourier.basis(rng,nbasis)
  }else if (type == "exp"){
    basis = create.exponential.basis(rng,nbasis)
  }else if (type == "const"){
    basis = create.constant.basis(rng)
  }else if (type == "mon"){
    basis = create.monomial.basis(rng,nbasis)
  }
  
  if(is.list(data)){
    data.list = data
  }else{
    data.list = sparse.list(data, time)
  }
  data.mat = do.call(rbind,data.list)
  
  if(!is.null(lambda)){
    curv.Lfd = int2Lfd(2)
    curv.fdPar = fdPar(basis,curv.Lfd,lambda)
    smooth = smooth.basis(data.mat[,1],data.mat[,-1],curv.fdPar)
  }else{
    smooth = smooth.basis(data.mat[,1],data.mat[,-1],basis)
  }
  
  if(ncol(smooth$fd$coefs)>1){
    smooth$fd$coefs = array(smooth$fd$coefs,dim=c(nrow(smooth$fd$coefs),1,ncol(smooth$fd$coefs)))
  }
  
  return(smooth$fd)
}

