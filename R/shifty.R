shifty <- function(x, y, shift)
{
#SHIFTY estimates value of Y for periodic data for
#       X shifted by amount SHIFT.
#  It is assumed that X spans interval over which functionis periodic.
#  Last modified 6 February 2001

ydim <- dim(y)
if (is.null(ydim)) ydim <- 1
if (length(ydim) > 3) stop("Y has more than three dimensions")

if (shift == 0) {
   yshift <- y
   return(yshift)
}

n   <- ydim[1]
xlo <- min(x)
xhi <- max(x)
wid <- xhi - xlo
if (shift > 0) {
   while (shift > xhi)  shift <- shift - wid 
   ind <- 2:n
   x2  <- c(x, x[ind]+wid)
   xshift <- x + shift
   if (length(ydim) == 1) {
	  y2 <- c(y, y[ind])
      yshift <- approx(x2, y2, xshift)$y
   } 
   if (length(ydim) == 2) {
	   nvar <- ydim[2]
	   yshift <- matrix(0,n,nvar)
      for (ivar in 1:nvar) {
         y2 <- c(y[,ivar], y[ind,ivar])
         yshift[,ivar] <- approx(x2, y2, xshift)$y
      }
   } 
   if (length(ydim) == 3) {
	   nrep <- ydim[2]
	   nvar <- ydim[3]
      yshift <- array(0,c(n,nrep,nvar))
      for (irep in 1:nrep) for (ivar in 1:nvar) {
         y2 <- c(y[,irep,ivar], y[ind,irep,ivar])
         yshift[,irep,ivar] <- approx(x2, y2, xshift)$y
      }
   } 
} else {
   while (shift < xlo - wid) shift <- shift + wid
   ind <- 1:(n-1)
   x2 <- c(x[ind]-wid, x)
   xshift <- x + shift
   if (length(ydim) == 1) {
      y2 <- c(y[ind], y)
      yshift <- approx(x2, y2, xshift)$y
   } 
   if (length(ydim) == 2) {
	   nvar <- ydim[2]
	   yshift <- matrix(0,n,nvar)
	   for (ivar in 1:nvar) {
		   y2 <- c(y[ind,ivar],y[,ivar])
		   yshift[,ivar] <- approx(x2, y2, xshift)$y
	   }
   }
   if (length(ydim) == 3) {
	   nrep <- ydim[2]
	   nvar <- ydim[3]
      yshift <- array(0, c(n,nrep,nvar))
      for (irep in 1:nrep) for (ivar in 1:nvar) {
         y2 <- c(y[ind,irep,ivar], y[,irep,ivar])
         yshift[,irep,ivar] <- approx(x2, y2, xshift)$y
      }
   }
}
return(yshift)
}
