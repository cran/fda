eigen.pda = function(pdaList,plotresult=TRUE,npts=501,...)
{
  rangval = pdaList$resfdlist[[1]]$basis$rangeval
  
  m = length(pdaList$resfdlist)  
  tfine = seq(rangval[1],rangval[2],length.out=npts)

  bwtlist = pdaList$bwtlist

  if(m == 1){
    d = length(bwtlist)
    xlabstr = names(bwtlist[[1]]$fd$fdnames)[[1]]

    betamat = array(0,c(npts,d,d))
      
    for(i in 1:d){
      betamat[,1,d-i+1] = -eval.fd(tfine,bwtlist[[i]]$fd)
      if(i < d) betamat[,i+1,i] = 1
    }
  }
  else{
    d = length(bwtlist[[1]][[1]])
    xlabstr = names(bwtlist[[1]][[1]][[1]]$fd$fdnames)[[1]]
#    betamat = array(0,c(npts,m,m,d))

    betamat = array(0,c(npts,m*d,m*d))
    
    for(k in 1:d){
      for(j in 1:m){
        for(i in 1:m){
#                betamat[,i,j,k] = eval.fd(tfine,bwtlist[[i]][[j]][[k]]$fd)
          betamat[,j,m*(d-k)+i] = -eval.fd(tfine,bwtlist[[i]][[j]][[k]]$fd)
        }
        if(k < d){
            betamat[,m*k+j,m*(k-1)+j] = 1
        }
      }
    }
  } 
  
  eigvals = matrix(0,npts,m*d)
              
  for(i in 1:npts){
      eigvals[i,] = eigen(betamat[i,,])$values
  }
  
  if(plotresult){
     par(mfrow=c(2,1))
     matplot(tfine,Re(eigvals),type='l',xlab=xlabstr,ylab='Real',main='Eigenvalues',...)
     abline(h = 0)
     matplot(tfine,Im(eigvals),type='l',xlab=xlabstr,ylab='Imaginary',...)
     abline(h = 0)
  }
  
  return(list(argvals=tfine,eigvals=eigvals))  
} 
  
  