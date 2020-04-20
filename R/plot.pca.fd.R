plot.pca.fd <- function(x, nx = 128, pointplot = TRUE, harm = 0,
                        expand = 0, cycle = FALSE, ...)
{
#
#	Plots the harmonics produced by PCA.FD.
#
#  If pointplot=TRUE, then the harmonics are plotted as + and -
#  otherwise lines are used.	Another thing that needs doing is an
#		 arrowplot option.
#
# If harm = 0 (the default) then all the computed harmonics are plotted.
#	 Otherwise those in jharm are plotted.
# If expand =0 then effect of +/- 2 standard deviations of each pc are given
#	 otherwise the factor expand is used.
# If cycle=TRUE and there are 2 variables then a cycle plot will be drawn
#	If the number of variables is anything else, cycle will be ignored.
#

# last modified 2007 May 3 by Spencer Graves	
#	previously modified 20 March 2006

  pcafd <- x
  if (!(inherits(pcafd, "pca.fd"))) 
    stop("Argument 'x' is not a pca.fd object.")

  harmfd	<- pcafd[[1]]
  basisfd <- harmfd$basis
  rangex	<- basisfd$rangeval
  {
    if(length(nx)>1){
      x <- nx
      nx <- length(x)
    }
    else    
      x <- seq(rangex[1], rangex[2], length = nx)
  }
  fdmat	 <- eval.fd(x, harmfd)
  meanmat <- eval.fd(x, pcafd$meanfd)
  dimfd	 <- dim(fdmat)
  nharm	 <- dimfd[2]
#
# check number of panels
  plotsPerPg <- sum(par("mfrow"))
#   
  harm	<- as.vector(harm)
  if(harm[1] == 0) harm <- (1:nharm)
#  
  if(length(dimfd) == 2) {
    for(jharm in 1:length(harm)) {
#    for(iharm in harm) {
      if(jharm==2){
        op <- par(ask=TRUE)
        on.exit(par(op)) 
      }
#        
      iharm <- harm[jharm] 
      if(expand == 0) {
        fac <- sqrt(pcafd$values[iharm])
      } else {
        fac <- expand
      }
#      
      vecharm <- fdmat[, iharm]
      pcmat <- cbind(meanmat + fac * vecharm, meanmat - fac * vecharm)
      if (pointplot) plottype <- "p" else plottype <- "l"
      percentvar <- round(100 * pcafd$varprop[iharm], 1)
      plot(x, meanmat, type = "l", ylim=c(min(pcmat),max(pcmat)),
           ylab=paste("Harmonic", iharm),
           main=paste("PCA function", iharm,
             "(Percentage of variability", percentvar, ")"),
           ...)
      if (pointplot) {
        points(x, pcmat[,1], pch="+")
        points(x, pcmat[,2], pch="-")
      } else {
        lines(x, pcmat[,1], lty=2)
        lines(x, pcmat[,2], lty=3)
      }
    }
  } else {
    if(cycle && dimfd[3] == 2) {
      meanmat <- drop(meanmat)
      for(jharm in 1:length(harm)) {
#      for(iharm in harm) {
        if(jharm==2){
          op <- par(ask=TRUE)
          on.exit(par(op)) 
        }
        iharm <- harm[jharm]
#
        {
          if(expand == 0) fac <- 2 * sqrt(pcafd$values[iharm])
          else fac <- expand
        }
        matharm <- fdmat[, iharm,	]
        mat1 <- meanmat + fac * matharm
        mat2 <- meanmat - fac * matharm
        if (pointplot) plottype <- "p" else plottype <- "l"
        percentvar <- round(100 * pcafd$varprop[iharm],1)
        plot(meanmat[,1], meanmat[,2], type=plottype,
             xlim=c(min(c(mat1[,1],mat2[,1])),max(c(mat1[,1],mat2[,1]))),
             ylim=c(min(c(mat1[,2],mat2[,2])),max(c(mat1[,2],mat2[,2]))),
             main=paste("PCA function", iharm,
               "(Percentage of variability", percentvar, ")"),
             ...)
        if (pointplot) {
          points(mat1[, 1], mat1[, 2], pch="+")
          points(mat2[, 1], mat2[, 2], pch="-")
        }
        else {
          lines (mat1[, 1], mat1[, 2], lty=2)
          lines (mat2[, 1], mat2[, 2], lty=3)
        }
      }
    }
    else {
      for(jharm in 1:length(harm)) {
#      for(iharm in harm) {
        if(jharm==2){
          op <- par(ask=TRUE)
          on.exit(par(op)) 
        }
        iharm <- harm[jharm]
#        
        fac <- {
          if (expand == 0) sqrt(pcafd$values[iharm]) 
          else expand
        }
        
        meanmat <- drop(meanmat)
        matharm <- fdmat[, iharm, ]
        nvar <- dim(matharm)[2]
        for (jvar in 1:nvar) {
          pcmat <- cbind(meanmat[, jvar] + fac * matharm[, jvar],
                         meanmat[, jvar] - fac * matharm[, jvar])
          if (pointplot) plottype <- "p" else plottype <- "l"
          percentvar <- round(100 * pcafd$varprop[iharm], 1)
          plot(x, meanmat[,jvar], type=plottype,
               ylab=paste("Harmonic", iharm),
               sub = paste("PCA function", iharm,
                 "(Percentage of variability", 
                 percentvar,")"),
               main = dimnames(fdmat)[[3]][jvar],
               ...)
          if (pointplot) {
            points(x, pcmat[,1], pch="+")
            points(x, pcmat[,2], pch="-")
          }
          else {
            lines (x, pcmat[,1], lty=2)
            lines (x, pcmat[,2], lty=3)
          }
        }
      }
    }
  }
  invisible(NULL)
}
