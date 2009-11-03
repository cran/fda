
#BD2
BD2=function(matrizDatos){
	n=dim(matrizDatos)[1]
	p=dim(matrizDatos)[2]
	cont=rep(0,n)
	for (i in 1:(n-1)){
		for (j in (i+1):n){
			cont=cont+estaEntre(c(i,j),matrizDatos)
		}
	}
	contg=(cont/combinat(n,2))
}

#indicator function
estaEntre=function(v,matrizDatos){
	n=dim(matrizDatos)[1]
	p=dim(matrizDatos)[2]
	Z=matrizDatos
	inf=t(apply(Z[v,],2,min))
	sup=t(apply(Z[v,],2,max))
	resultados=colSums((t(Z)<=t(sup)%*%rep(1,n))* (t(Z)>=t(inf)%*%rep(1,n)))==p
}

#combination
combinat=function(n,p){
	if (n<p){combinat=0}
	else {combinat=factorial(n)/(factorial(p)*factorial(n-p))}
}

#BD3
BD3=function(matrizDatos){
	n=dim(matrizDatos)[1]
	p=dim(matrizDatos)[2]
	cont=rep(0,n)
	for (i in 1:(n-2)){
		for (j in (i+1):(n-1)){
			for (k in (j+1):n){
				cont=cont+estaEntre(c(i,j,k),matrizDatos)
			}
		}
	}
	contg=(cont/combinat(n,3))
}

#MBD
MBD=function(matrizDatos){
	n=dim(matrizDatos)[1]
	p=dim(matrizDatos)[2]
	cont=rep(0,n)
	for (i in 1:(n-1)){
		for (j in (i+1):n){
			cont=cont+aprops(c(i,j),matrizDatos)
		}
	}
	contg=(cont/combinat(n,2))
}

#proportion function
aprops=function(v,matrizDatos){
	n=dim(matrizDatos)[1]
	p=dim(matrizDatos)[2]
  Z=matrizDatos
	inf=t(apply(Z[v,],2,min))
	sup=t(apply(Z[v,],2,max))
	resul=colSums((t(Z)<=t(sup)%*%rep(1,n))* (t(Z)>=t(inf)%*%rep(1,n)))
	resultado=(resul/p)
}




#function boxplot
#fit: p by n functional data matrix, n is the number of curves
#method: BD2, BD3, MBD
fbplot=function(fit,x=NULL,method='MBD',depth=NULL,plot=TRUE,prob=0.5,color=6,outliercol=2,
				barcol=4,fullout=FALSE, factor=1.5,...){
				
  if(is.fdSmooth(fit) | is.fdPar(fit)){ fit = fit$fd }  
	if(is.fd(fit)){
    if(length(x)==0){
      x = seq(fit$basis$rangeval[1],fit$basis$rangeval[2],len=101)
    }
    fit = eval.fd(x,fit)
  }				
				
	tp=dim(fit)[1]
	n=dim(fit)[2]
	if (length(x)==0) {x=1:tp}
  #compute band depth	
  if (length(depth)==0){
	if (method=='BD2') {depth=BD2(t(fit))}
	else if (method=='BD3') {depth=BD3(t(fit))}
	else if (method=='MBD') {depth=MBD(t(fit))}
	else if (method=='Both') {depth=round(BD2(t(fit)),4)*10000+MBD(t(fit))}
  }

	dp_s=sort(depth,decreasing=T)
	index=order(depth,decreasing=T)
	if (plot) {
	plot(x,fit[,index[1]],lty=1,lwd=2,col=1,type='l',...)
	}
	for (pp in 1:length(prob)){
		m=ceiling(n*prob[pp])#at least 50%
		center=fit[,index[1:m]]
		out=fit[,index[(m+1):n]]
		inf=apply(center,1,min)
		sup=apply(center,1,max)
		
		if (prob[pp]==0.5){ #check outliers
			dist=factor*(sup-inf)
			upper=sup+dist
			lower=inf-dist
			outly=(fit<=lower)+(fit>=upper)
			outcol=colSums(outly)
			remove=(outcol>0)
			#outlier column
			colum=1:n
			outpoint=colum[remove==1]
			out=fit[,remove]
			woout=fit
			good=woout[,(remove==0),drop=FALSE]
			maxcurve=apply(good,1,max)
			mincurve=apply(good,1,min)
			if (sum(outly)>0){
				if (plot) {
				matlines(x,out,lty=2,col=outliercol,type='l',...)
				}
			}
			barval=(x[1]+x[tp])/2
			bar=which(sort(c(x,barval))==barval)[1]
			if (plot) {
			lines(c(x[bar],x[bar]),c(maxcurve[bar],sup[bar]),col=barcol,lwd=2)
		    lines(c(x[bar],x[bar]),c(mincurve[bar],inf[bar]),col=barcol,lwd=2)
			}
		}
		xx=c(x,x[order(x,decreasing=T)])
		supinv=sup[order(x,decreasing=T)]
		yy=c(inf,supinv)
		if (plot) {
		if (prob[pp]==0.5) {polygon(xx,yy,col=color[pp],border=barcol,lwd=2)}
		else {polygon(xx,yy,col=color[pp],border=NA)}
		}
	}
	if (plot) {
	lines(x,fit[,index[1]],lty=1,lwd=2,col=1,type='l')
	lines(x,maxcurve,col=barcol,lwd=2)
	lines(x,mincurve,col=barcol,lwd=2)
	if (fullout) {
		if (sum(outly)>0){
				if (plot) {
				matlines(x,out,lty=2,col=outliercol,type='l',...)
				}
			}
		}
	}
	return(list(depth=depth,outpoint=outpoint))
}


