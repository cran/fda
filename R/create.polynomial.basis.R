create.polynomial.basis <- function (rangeval=c(0,1), nbasis=2, ctr=0,
                   dropind=NULL, quadvals=NULL, values=NULL,
                   basisvalues=NULL, names='polynom', axes=NULL)
{

#  This function creates a polynomial functional data basis, for
#    polynomials of the form  (x - c)^j
#  Arguments
#  RANGEVAL ... an array of length 2 containing the lower and upper
#               boundaries for the rangeval of argument values
#  NBASIS   ... the number of basis functions
#  CTR      ... The centering constant C.  By default, this is 0.
#  DROPIND  ... A vector of integers specificying the basis functions to
#               be dropped, if any.
#  Returns
#  BASISOBJ ... a functional data basis object of type "polynomial"

#  Last modified 9 November 2008 by Spencer Graves
#  Last modified 3 January 2008 by Jim Ramsay

##
## 1.  check RANGEVAL
##
  if(!is.numeric(rangeval))
    stop('rangaval must be numeric;  class(rangeval) = ',
         class(rangeval) )
  if(length(rangeval)<1)
    stop('rangeval must be a numeric vector of length 2;  ',
         'length(rangeval) = 0.')
  if (length(rangeval) == 1) {
    if( rangeval <= 0)
      stop("rangeval a single value that is not positive:  ",
           rangeval)
    rangeval <- c(0,rangeval)
  }
  if(length(rangeval)>2)
    stop('rangeval must be a vector of length 2;  ',
         'length(rangeval) = ', length(rangeval))
  if(diff(rangeval)<=0)
    stop('rangeval must cover a positive range;  diff(rangeval) = ',
         diff(rangeval) )
##
## 2.  check nbasis>0
##
#  check NBASIS
  if(!is.numeric(nbasis))
    stop('nbasis must be numeric;  class(nbasis) = ',
         class(nbasis) )
  if(length(nbasis) != 1)
    stop('nbasis must be a scalar;  length(nbasis) = ',
         length(nbasis) )
  if((nbasis%%1) != 0)
    stop('nbasis must be an integer;  nbasis%%1 = ',
         nbasis%%1)
  nbasis <- as.integer(nbasis)
  if (nbasis <= 0) stop ("NBASIS is not positive.")
##
## 3.  check DROPIND
##
  if (length(dropind) == 0) dropind <- NULL
#
  if (length(dropind) > 0){
    if(!is.numeric(dropind))
      stop('dropind must be numeric;  is ', class(dropind))
    doops <- which((dropind%%1)>0)
    if(length(doops)>0)
      stop('dropind must be integer;  element ', doops[1],
           " = ", dropind[doops[1]], '; fractional part = ',
           dropind[doops[1]] %%1)
#
    doops0 <- which(dropind<=0)
    if(length(doops0)>0)
      stop('dropind must be positive integers;  element ',
           doops0[1], ' = ', dropind[doops0[1]], ' is not.')
    doops2 <- which(dropind>nbasis)
    if(length(doops2)>0)
        stop("dropind must not exceed nbasis = ", nbasis,
             ';  dropind[', doops2[1], '] = ', dropind[doops2[1]])
#
    dropind <- sort(dropind)
    if(length(dropind) > 1) {
      if(min(diff(dropind)) == 0)
        stop("Multiple index values in DROPIND.")
    }
  }
##
## 5.  Check ctr
##
  if(!is.numeric(ctr))
    stop('ctr must be numeric;  class(ctr) = ',
         class(ctr) )
  if(length(ctr) != 1)
    stop('ctr must be a scalar;  length(ctr) = ',
         length(ctr) )
##
## 6.  set up the basis object
##
  type        <- "polynom"
  params      <- as.vector(ctr)
#
  basisobj <- basisfd(type=type,   rangeval=rangeval, nbasis=nbasis,
                    params=params, dropind=dropind,   quadvals=quadvals,
                    values=values, basisvalues=basisvalues)
##
## 7.  names
##
  {
    if(length(names) == nbasis)
      basisobj$names <- names
    else {
      if(length(names)>1)
        stop('length(names) = ', length(names), ';  must be either ',
             '1 or nbasis = ', nbasis)
      basisobj$names <- paste(names, 0:(nbasis-1), sep="")
    }
  }
##
## 8.  Done
##
  if(!is.null(axes))basisobj$axes <- axes
  basisobj

}
