#Fperm.fd(annualprec, xfdlist, betalist, nperm=10)

Fperm.fd <- function(yfdPar, xfdlist, betalist, wt=NULL,
            nperm=200,argvals=NULL,q=0.05,plotres=TRUE,...)
{
# yfdpar, xfdList, betalist, wt = standard inputs to fRegress
# nperm = number of permutations,
# argvals = where to evaluate functional responses,
# q =  quantile to compare
# plotres:  Do we plot the results?

  # set up containers for null results
  
  Fnull     = rep(0,nperm)
  Fnullvals = c()
  
  # switch q from residual to main body
  
  q <- 1-q
  
  # keep track of processing time
  
  begin <- proc.time()
  fRegressList <- fRegress(yfdPar, xfdlist, betalist)
  elapsed.time <- max(proc.time()-begin,na.rm=TRUE)
  if( elapsed.time > 30/nperm ) {
    print(paste('Estimated Computing time =',
                round(nperm*elapsed.time),'seconds.'))
  }
  
  #. obtain fit vector yhat from fRegressList
  
  yhat <- fRegressList$yhatfdobj. #. the fitted values
  if(is.null(yhat)) yhat = fRegressList$yhat
  if(is.list(yhat) && ('fd' %in% names(yhat))) yhat <- yhat$fd
  
  #. initial invokation of Fstat.fd
  
  tFstat <- Fstat.fd(yfdPar,yhat,argvals)
  
  # two sets of results 
  
  Fvals <- tFstat$F
  Fobs  <- max(Fvals)
  
  #.  obtain argument values
  
  argvals <- tFstat$argvals
  
  #. determine length of data vector
  
  if(is.vector(yfdPar)){ 
    n <- length(yfdPar) 
  }  else { 
    n <- ncol(yfdPar$coefs) 
  }
  
  #. --------------------------------------------------------------------------
  #                   Loop through permutations
  #. --------------------------------------------------------------------------
  
  for(i in 1:nperm){
    
    tyfdPar <- yfdPar[sample(n)]
    
    #  analyze of a vector of observations 
    
    fRegressList <- fRegress(tyfdPar, xfdlist, betalist)
    
    #  extract fit vector from fRegressList
    
    if(is.fd(yhat)) {
      yhat <- fRegressList$yhatfdobj
      if(is.list(yhat) && ('fd' %in% names(yhat))) yhat <- yhat$fd
    } else { 
      yhat <- fRegressList$yhat 
    }
    
    # ith analysis by Fstat.fd
    
    tFstat <- Fstat.fd(tyfdPar,yhat,argvals)
    
    # add Fnullvals to container Fnullvals
    
    Fnullvals <- cbind(Fnullvals,tFstat$F)
    
    #. 
    Fnull[i] <- max(Fnullvals[,i])
  }
  
  #. --------------------------------------------------------------------------
  #                   Display results
  #. --------------------------------------------------------------------------

  pval <- mean( Fobs < Fnull )
  qval <- quantile(Fnull,q)
  
  pvals.pts <- apply(Fvals<Fnullvals,1,mean)
  qvals.pts <- apply(Fnullvals,1,quantile,q)
  
  if(plotres) {
    #. dispslay results using scattter plots
    if(is.fd(yfdPar)){
      ylims <- c(min(c(Fvals,qval,qvals.pts)),max(c(Fobs,qval)))
      
      if( is.null(names(yhat$fdnames)) ){ xlab <- 'argvals' }
      else{ xlab <- names(yhat$fdnames)[1] }
      
      plot(argvals,Fvals,type="l",ylim=ylims,col=2,lwd=2,
           xlab=xlab,ylab='F-statistic',main='Permutation F-Test',...)
      lines(argvals,qvals.pts,lty=3,col=4,lwd=2)
      abline(h=qval,lty=2,col=4,lwd=2)
      
      legendstr <- c('Observed Statistic',
                    paste('pointwise',1-q,'critical value'),
                    paste('maximum',1-q,'critical value'))
      
      legend(argvals[1],ylims[2],legend=legendstr,col=c(2,4,4),
             lty=c(1,3,2),lwd=c(2,2,2))
    } else {
      #. display results with histogram
      
      xlims <- c(min(c(Fnull,Fobs)),max(c(Fnull,Fobs)))
      hstat <- hist(Fnull,xlim=xlims,lwd=2,xlab='F-value',
                   main = 'Permutation F-Test')
      abline(v = Fobs,col=2,lwd=2)
      abline(v = qval,col=4,lty=2,lwd=2)
      
      legendstr <- c('Observed Statistic',
                    paste('Permutation',1-q,'critical value'))
      
      legend(xlims[1],max(hstat$counts),legend=legendstr,col=c(2,4),
             lty=c(1,2),lwd=c(2,2))
    }
  }
  
  #. return results
  
  return(list(pval=pval, qval=qval, Fobs=Fobs, Fnull=Fnull,
              Fvals=Fvals, Fnullvals=Fnullvals, pvals.pts=pvals.pts,
              qvals.pts=qvals.pts,
              fRegressList=fRegressList, argvals=argvals))
}
