sparse.list <- function(data,time){
  #Create a list with sparse data from a matrix that has NA
  #
  # Arguments:
  #
  # DATA .... Sparse data -- If the set is supplied as a matrix object, the rows must correspond to argument values and 
  #           columns to replications, and it will be assumed that there is only one variable per observation. If y is a 
  #           three-dimensional array, the first dimension corresponds to argument values, the second to replications, 
  #           and the third to variables within replications. 
  # TIME .... Time points where the observations where taken. length(time) == nrow(data)
  
  ndim = length(dim(data))
  if(ndim ==3){
    ind = apply(data[,,1],2, function(x) which(!(is.na(x))))
    y = lapply(1:dim(data)[2], function(x) cbind(time[ind[[x]]],data[ind[[x]],x,]))
  }else{
    ind = apply(data,2, function(x) which(!(is.na(x))))
    y = lapply(1:dim(data)[2], function(x) cbind(time[ind[[x]]],data[ind[[x]],x]))
  }
  
  
  return(y)
}
