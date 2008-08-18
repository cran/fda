as.array3 <- function(x){
  dimx <- dim(x)
  ndim <- length(dimx)
#   
  if(ndim==3)return(x)
#  
  xName <- substring(deparse(substitute(x)), 1, 22) 
  if(ndim>3)
    stop('length(dim(', xName, ") = ", ndim,
         ' > 3')
#
  x. <- as.matrix(x)
  xNames <- dimnames(x.)
  dim(x.) <- c(dim(x.), 1)
#
  if(is.list(xNames))
    dimnames(x.) <- list(xNames[[1]], xNames[[2]], NULL)
#
  x. 
}

