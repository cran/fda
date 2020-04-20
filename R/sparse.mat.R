sparse.mat <- function(datalist){
  #
  #Create a matrix of sparse data with NAs out of a list of sparse data
  #
  # Arguments:
  #
  # DATALIST .... A list object. Each element of the list is a matrix with ncol > 1. 
  #               The first column of each element corresponds to the point index per 
  #               observation. The following columns are the observations per variable.
  
  time = sort(unique(unlist(lapply(datalist, function(x) x[,1]))))
  data = matrix(NA, nrow = length(time), ncol = length(datalist))
  
  nvar = ncol(datalist[[1]]) - 1
  
  if(nvar == 1){
    data = matrix(NA, nrow = length(time), ncol = length(datalist))
    for(i in 1:length(datalist)){
      data[which((time %in% datalist[[i]][,1])),i] = datalist[[i]][,2]
    }
  }else{
    data = array(NA, dim = c(length(time),length(datalist),nvar))
    for(j in 1:nvar){
      for(i in 1:length(datalist)){
        data[which((time %in% datalist[[i]][,1])),i,j] = datalist[[i]][,2]
      }
    }
  }
  return(data)
}
