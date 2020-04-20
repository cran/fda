hex <- function(ctr=c(0,0), rad=1, sig=0, meshlevel=1) {
  #  HEX sets up six boundary points around a hexagon,
  #  the PET architecture, and data at 19 locations
  angle <- seq(0,2*pi,len=7)
  x <- round(rad*cos(angle) + ctr[1],7)
  y <- round(rad*sin(angle) + ctr[2],7)
  edg <- matrix(0,7,2)
  edg[1:6,] <- cbind(x[1:6],y[1:6])
  edg[  7,] <- edg[1,]
  #  define the mesh
  if (meshlevel == 1) {
    pts  <- matrix(0,7,2)
    tri  <- matrix(0,6,3)
    loc  <- matrix(0,19,2)
    pts[1:6,] <- edg[1:6,]
    tri[,1] <- 7
    tri[1:6,2]  <- 1:6
    tri[   ,3]  <- c(2:6,1)
    loc[ 1: 6,] <- pts[1:6,]
    loc[ 7:12,] <- (pts[1:6,1] + pts[1:6,2])/2
    loc[13:19,] <- pts/2
    edg <- edg[1:6,]
  } else if (meshlevel == 2) {
    pts <- matrix(0,19,2)
    #  set up points matrix
    mpts2 <- 0
    for (i in 1:6) {
      mpts1 <- mpts2 + 1
      mpts2 <- mpts2 + 3
      pts[mpts1,  ] <-  edg[i,]
      pts[mpts1+1,] <-  edg[i,]/2
      pts[mpts1+2,] <- (edg[i,] + edg[i+1,])/2
    }
    pts[   19,] <- c(0,0)
    #  set up triangle matrix
    tri <- matrix(0,24,3)
    mpts2 <- 0
    mtri2 <- 0
    for (i in 1:6) {
      mpts1 <- mpts2 + 1
      mpts2 <- mpts2 + 3
      mtri1 <- mtri2 + 1
      mtri2 <- mtri2 + 4
      
      tri[mtri1  ,1] <- mpts1
      tri[mtri1  ,2] <- mpts1+2
      tri[mtri1  ,3] <- mpts1+1
      
      tri[mtri1+1,1] <- mpts1+1
      tri[mtri1+1,2] <- mpts1+2
      tri[mtri1+1,3] <- (mpts1+4) %% 18
      
      tri[mtri1+2,1] <- mpts1+2
      tri[mtri1+2,2] <- (mpts1+3) %% 18
      tri[mtri1+2,3] <- (mpts1+4) %% 18
      
      tri[mtri1+3,1] <- mpts1+1
      tri[mtri1+3,2] <- (mpts1+4) %% 18
      tri[mtri1+3,3] <- 19        
    }
    pts <- rbind(pts[1:18,], matrix(0,1,2))
    edg <- edg[1:6,]
    loc <- pts
  } else {
    stop('Argument meshlevel is neither 1 or 2.')
  }
  dat <- 1 - loc[,1]^2 - loc[,2]^2
  dat <- dat + rnorm(19,1)*sig
  return(list(p=pts, e=edg, t=tri, loc=loc, dat=dat))
}