scoresPACE <- function(data, time, covestimate, PC){
  if(is.list(data)){
    datamat = sparse.mat(data)
    datalist = sparse.list(datamat, time)
  }else{
    datamat = data
    datalist = sparse.list(datamat, time)
  }
  

  indexes = lapply(datalist, function(x) which(time %in% x[,1])) #sampling points for each subject
  coef = PC$values
  phi = lapply(indexes, function(x) eval.fd(x,PC$harmonics))
  
  
  if(class(covestimate$cov.estimate) == "bifd"){
    nvar = 1
    mean.point = matrix(eval.fd(time, covestimate$meanfd),nrow = nrow(eval.fd(time, covestimate$meanfd)),ncol = nvar)
    data.var.mat = (apply(datamat,2, function(x) x - mean.point))^2 #data given as matrix
    mu = lapply(indexes,function(x) as.matrix(eval.fd(x,covestimate$meanfd),nrow = length(x),ncol = nvar))
    cent = lapply(1:length(datalist), function(x) datalist[[x]][,-1]-mu[[x]])
    
    covdiag = diag(eval.bifd(time,time,covestimate$cov.estimate))
    varest.mat = as.matrix(apply(data.var.mat,2, function(x) x-covdiag))
    varest = mean(as.vector(varest.mat)[!is.na(as.vector(varest.mat))])
    diagon = lapply(1:length(datalist), function(x) diag(varest, nrow = length(indexes[[x]]), ncol = length(indexes[[x]])))
    sigmai = Map('+', lapply(indexes, function(x) eval.bifd(x,x,covestimate$cov.estimate)), diagon)

    scr = t(do.call(cbind,lapply(1:length(datalist), function(x) coef*(t(phi[[x]])%*%solve(sigmai[[x]])%*%cent[[x]]) )))
    
  }else{
    k=0 #nvar
    l = length(covestimate$cov.estimate)
    while(l>0){
      k = k + 1
      l = l - k
    }
    nvar = k 
    mean.point = matrix(eval.fd(time, covestimate$meanfd),nrow = nrow(eval.fd(time, covestimate$meanfd)),ncol = nvar)
    data.var.mat = (apply(datamat,2, function(x) x - mean.point))^2 #data given as matrix
    mu = lapply(indexes,function(x) matrix(eval.fd(x,covestimate$meanfd),nrow = length(x),ncol = nvar))
    cent = lapply(1:length(datalist), function(x) datalist[[x]][,-1]-mu[[x]])
    
    ind = cumsum(c(1,k - 0:(k-2)))
    covdiag = as.vector(do.call(cbind, lapply(covestimate$cov.estimate[ind], function(x) diag(eval.bifd(time,time,x)))))
    varest.mat = as.matrix(apply(data.var.mat,2, function(x) x-covdiag))
    ntime = dim(datamat)[1]
    varest = unlist(lapply(1:nvar, function(x) mean(as.vector(varest.mat[(ntime*x-(ntime-1)):(ntime*x),])[!is.na(as.vector(varest.mat[(ntime*x-(ntime-1)):(ntime*x),]))])))
    
    diagon = list()
    sigmai = list()
    scr = array(NA, dim = c(ncol(datamat),ncol(coef),nvar))
    for(l in 1:nvar){
      diagon[[l]] = lapply(1:ncol(datamat), function(x) diag(varest[l], nrow = length(indexes[[x]]), ncol = length(indexes[[x]])))
      sigmai[[l]] = Map('+', lapply(indexes, function(x) eval.bifd(x,x,covestimate$cov.estimate[[ind[l]]])), diagon[[l]])
      scr[,,l] = t(do.call(cbind,lapply(1:ncol(datamat), function(x) coef[l,]*(t(phi[[x]][,,l])%*%solve(sigmai[[l]][[x]])%*%cent[[x]][,l]) )))
    }
  }
  
  
  return(scr)
}
