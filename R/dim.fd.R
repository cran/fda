dim.fd <- function(x)
{
  if(!is.fd(x))
    return(dim(x)) 
  else {
    d = dim(x$coefs)
    return( d[2:length(d)] )
  }
}
