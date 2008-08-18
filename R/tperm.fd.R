tperm.fd <- function(x1fd,x2fd,nperm=200,q=0.95,argvals=NULL,plotres=TRUE) # first and second 
{                                                                          # groups of data,
    if( !is.fd(x1fd) | !is.fd(x2fd) ){                                     # number permuts
        stop("x1fd and x2fd must both be functional data objects")         # quantile
    }                                                                      # where to evaluate
                                                                           # do I plot
    rangeobs = x1fd$basis$range
    rangehat = x2fd$basis$range

    if( !prod(rangeobs == rangehat) ){
        stop("x1fd and x2fd do not have the same range.")
    }

    if(is.null(argvals)){
        argvals = seq(rangeobs[1],rangeobs[2],length.out=101)
    }

    x1mat = eval.fd(argvals,x1fd)
    x2mat = eval.fd(argvals,x2fd)

    n1 = ncol(x1mat)
    n2 = ncol(x2mat)

    Xmat = cbind(x1mat,x2mat)

    Tnull = rep(0,nperm)

    Tnullvals = matrix(0,length(argvals),nperm)

    for(i in 1:nperm){
        tXmat = Xmat[,sample(n1+n2)]

        tmean1 = apply(tXmat[,1:n1],1,mean)
        tmean2 = apply(tXmat[,n1+(1:n2)],1,mean)

        tvar1 = apply(tXmat[,1:n1],1,var)/n1
        tvar2 = apply(tXmat[,n1+(1:n2)],1,var)/n2

        Tnullvals[,i] = abs(tmean1-tmean2)/sqrt(tvar1+tvar2)
        Tnull[i] = max(Tnullvals[,i])
    }

    mean1 = apply(Xmat[,1:n1],1,mean)
    mean2 = apply(Xmat[,n1+(1:n2)],1,mean)

    var1 = apply(Xmat[,1:n1],1,var)/n1
    var2 = apply(Xmat[,n1+(1:n2)],1,var)/n2

    Tvals = abs(mean1-mean2)/sqrt(var1+var2)
    Tobs = max(Tvals)

    pval = mean( Tobs < Tnull )
    qval = quantile(Tnull,q)

    pvals.pts = apply(Tvals<Tnullvals,1,mean)
    qvals.pts = apply(Tnullvals,1,quantile,q)

    if(plotres){
        ylims = c( min(Tobs,-qval),max(Tobs,qval))

        plot(argvals,Tvals,type='l',col=2,ylim=ylims)
        lines(argvals,qvals.pts,lty=2,col=4)
        abline(h=qval,lty=2,col=4)
        abline(h=-qval,lty=2,col=4)
        abline(h=Tobs,lty=2,col=2)
    }


    return( list(pval=pval,qval=qval,Tobs=Tobs,Tnull=Tnull,
        Tvals=Tvals,Tnullvals=Tnullvals,qvals.pts=qvals.pts,
        pvals.pts=pvals.pts,argvals=argvals) )
}
