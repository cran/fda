triplot <- function(p, t) {
  #  p ... NP by 2 matrix of point coordinates
  #  t ... NT by 3 (or 4) matrix of indices of triangle vertices in P
  
  #  Last modified 12 June 2015 by Jim Ramsay
  
  nt = dim(t)[1]
  np = dim(p)[1]
  plot(p, lwd=2, xlab="w", ylab="s")
  for (i in 1:nt) {
    lines(c(p[t[i,1],1],p[t[i,2],1]), c(p[t[i,1],2],p[t[i,2],2]), lwd=2)
    lines(c(p[t[i,2],1],p[t[i,3],1]), c(p[t[i,2],2],p[t[i,3],2]), lwd=2)   
    lines(c(p[t[i,3],1],p[t[i,1],1]), c(p[t[i,3],2],p[t[i,1],2]), lwd=2)    
  } 
}
