wtcheck = function(n, wtvec) {

  if (is.null(wtvec)) return(rep(1,n))

  if (length(wtvec) != n) stop("'wtvec' of wrong length")

  if(!is.numeric(wtvec))
    stop("'wtvec' must be numeric;  class(wtvec) = ", class(wtvec))
  wdim <- dim(wtvec)
  if(length(wdim)>0){
    if(wdim[1] != n)
      stop("'wtvec' must be a column vector;  dim(wtvec) = ",
           paste(wdim, collapse=', '))
    wtNms <- dimnames(wtvec)[[1]]
    if(is.null(wtNms)) wtNms <- names(wtvec)
    wtvec <- as.numeric(wtvec)
    names(wtvec) <- wtNms
  }

  if (min(wtvec) <= 0)stop("All values of 'wtvec' must be positive.")

  return(wtvec)

}
