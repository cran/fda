Fperm.fd <- function(yfdPar, xfdlist, betalist,wt=NULL, # Standard inputs to fRegress
            nperm=200,argvals=NULL,q=0.95,plotres=TRUE) # number of permutations,
{                                                       # where to evaluate functional
    Fnull = rep(0,nperm)                                # responses, quantile to compare
                                                        # and do we plot the results?
    Fnullvals = c()

    begin <- proc.time()
    fRegressList <- fRegress(yfdPar, xfdlist, betalist)
    elapsed.time <- max(proc.time()-begin,na.rm=TRUE)       

    if( elapsed.time > 30/nperm ){
        print(paste('Estimated Computing time =',round(nperm*elapsed.time),'seconds.'))
    }

    yhat <- fRegressList$yhatfdobj
        
    tFstat <- Fstat.fd(yfdPar,yhat,argvals)

    Fvals <- tFstat$F
    Fobs = max(Fvals)

    argvals = tFstat$argvals

    if(is.vector(yfdPar)){ n = length(yfdPar) }
    else{ n = ncol(yfdPar$coefs) }

    for(i in 1:nperm){
        
        tyfdPar = yfdPar[sample(n)]

        fRegressList <- fRegress(tyfdPar, xfdlist, betalist)
        
        yhat <- fRegressList$yhatfdobj

        tFstat = Fstat.fd(yfdPar,yhat,argvals)
        
        Fnullvals <- cbind(Fnullvals,tFstat$F)
        
        Fnull[i] = max(Fnullvals[,i])
    }
    

    pval = mean( Fobs < Fnull )
    qval = quantile(Fnull,q)

    pvals.pts = apply(Fvals<Fnullvals,1,mean)
    qvals.pts = apply(Fnullvals,1,quantile,q)


    if(plotres){
        if(is.fd(yfdPar)){
            ylims = c(min(c(Fvals,qval)),max(c(Fobs,qval)))
    
            plot(argvals,Fvals,type="l",ylim=ylims,col=2)
            lines(argvals,qvals.pts,lty=2,col=4)
            abline(h=qval,lty=2,col=4)
            abline(h=Fobs,col=2)
        }
        else{
            xlims = c(min(c(Fnull,Fobs)),max(c(Fnull,Fobs)))
            hist(Fnull,xlim=xlims)
            abline(v = Fobs,col=2)
            abline(v = qval,col=4)
        }
    }
    
    return(list(pval=pval,qval=qval,Fobs=Fobs,Fnull=Fnull,
        Fvals=Fvals,Fnullvals=Fnullvals,pvals.pts=pvals.pts,qvals.pts=qvals.pts,
        fRegressList=fRegressList,argvals=argvals))
}
