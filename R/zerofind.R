zerofind <- function(fmat)
{
  frng <- range(fmat)
  if (frng[1] <= 0 && frng[2] >= 0) zeroin <- TRUE else zeroin <- FALSE
  return(zeroin)
}

