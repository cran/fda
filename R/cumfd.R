cumfd <- function(xrnd, xrng, nbreaks=7, nfine=101) {
  #  Compute cdf_fd over a closed interval using smooth.morph.  
  #  Only the values of x within the interior of xrng are used 
  #  in order to avoid distortion due to boundary inflation.
  #  Arguments:
  #  xrnd    ... A vector of variable values (unsorted)
  #  xrng    ... A vector of length 2 containing the boundary values.
  #  Wnbasis ... Number of basis functions used by smooth.morph.
  #  Snbasis ... Number of basis functions used by smooth.basis.
  #  nfine   ... Length of vector of equally spaced values spanning xrng.
  
  #  Last modified 25 March 2022 by Jim Ramsay
  
  #  check that values of x are within xrng
  
  if (min(xrnd) < xrng[1] || max(xrnd) > xrng[2]) 
    stop("Values of x outside of xrng.")
  
  #  sort the data and set up probability values
  
  xsort  <- sort(xrnd[xrnd > xrng[1] & xrnd < xrng[2]])
  N      <- length(xsort)
  prbvec <- (1:N)/(N+1)
  
  #  add boundary values
  
  pmesh <- c(0, prbvec, 1)
  xmesh <- c(xrng[1], xsort,  xrng[2])
  
  #  set up fdPar object for smooth.morph
  
  index = c(1, round(N*2:(nbreaks-1)/nbreaks), N+2)
  
  Wnorder <- 4
  Wnbasis <- nbreaks + Wnorder - 2
  Wbreaks <- xmesh[index]
  Wbasis  <- create.bspline.basis(xrng, Wnbasis, Wnorder, Wbreaks)
  WfdPar  <- fdPar(fd(matrix(0,Wnbasis,1), Wbasis))
  
  #  use smooth.morph to map sorted data into the interior of [0,1]
  
  result  <- smooth.morph(xmesh, pmesh, c(0,1), WfdPar)
  xfine   <- seq(0,1,len=101)
  Wfdobj  <- result$Wfdobj
  
  cdffine <- result$hfine
  cdffine[1]               <- 0
  cdffine[length(cdffine)] <- 1
  
  # plot(xfine, cdffine, type="l")
  # points(xmesh, pmesh)
  
  # plot(Wfdobj)
  
  return(list(Wfdobj=Wfdobj, cdffine=cdffine))
  
}