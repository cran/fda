
#  Generator function of class basisfd

basisfd <- function(type, rangeval, nbasis, params, dropind=NULL, 
                  quadvals=NULL, values=vector("list",0))
{
#  BASISFD  generator function of "basisfd" class.
#  Arguments:
#  BASISTYPE ... a string indicating the type of basisobj$  This may be one of
#               "Fourier", "fourier", "Fou", "fou",
#               "Bspline", "bspline", "Bsp", "bsp",
#               "pol", "poly", "polynomial",
#               "mon", "monom", "monomial",
#               "con", "const", "constant"
#               "exp", "exponen", "exponential"
#               "polyg" "polygon", "polygonal"
#  RANGEVAL ... an array of length 2 containing the lower and upper
#               boundaries for (the rangeval of argument values
#  NBASIS   ... the number of basis functions
#  PARAMS   ... If the basis is "fourier", this is a single number indicating
#                 the period.  That is, the basis functions are periodic on
#                 the interval (0,PARAMS) or any translation of it.
#               If the basis is "bspline", the values are interior points at
#                 which the piecewise polynomials join.
#                 Note that the number of basis functions NBASIS is equal
#                 to the order of the Bspline functions plus the number of
#                 interior knots, that is the length of PARAMS.
#               This means that NBASIS must be at least 1 larger than the
#                 length of PARAMS.
#  DROPIND ... A set of indices in 1:NBASIS of basis functions to drop
#                when basis objects are arguments.  Default is NULL  Note
#                that argument NBASIS is reduced by the number of indices,
#                and the derivative matrices in VALUES are also clipped.
#  QUADVALS .. A NQUAD by 2 matrix.  The firs t column contains quadrature
#                points to be used in a fixed point quadrature.  The second
#                contains quadrature weights.  For example, for (Simpson"s
#                rule for (NQUAD = 7, the points are equally spaced and the
#                weights are delta.*[1, 4, 2, 4, 2, 4, 1]/3.  DELTA is the
#                spacing between quadrature points.  The default is NULL.
#  VALUES  ... A list, with entries containing the values of
#                the basis function derivatives starting with 0 and
#                going up to the highest derivative needed.  The values
#                correspond to quadrature points in QUADVALS and it is
#                up to the user to decide whether or not to multiply
#                the derivative values by the square roots of the
#                quadrature weights so as to make numerical integration
#                a simple matrix multiplication.
#                Values are checked against QUADVALS to ensure the correct
#                number of rows, and against NBASIS to ensure the correct
#                number of columns.
#                The default is VALUES{1} = NULL
#  Returns
#  BASISOBJ  ... a basisfd object with slots
#         type
#         rangeval
#         nbasis
#         params
#         dropind
#         quadvals
#         values
#  Slot VALUES contains values of basis functions and derivatives at
#   quadrature points weighted by square root of quadrature weights.
#   These values are only generated as required, and only if slot
#   quadvals is not empty.
#
#  An alternative name for (this function is CREATE_BASIS, but PARAMS argument
#     must be supplied.
#  Specific types of bases may be set up more conveniently using functions
#  CREATE_BSPLINE_BASIS     ...  creates a b-spline basis
#  CREATE_CONSTANT_BASIS    ...  creates a constant basis
#  CREATE_EXPONENTIAL_BASIS ...  creates an exponential basis
#  CREATE_FOURIER_BASIS     ...  creates a fourier basis
#  CREATE_MONOMIAL_BASIS    ...  creates a monomial basis
#  CREATE_POLYGON_BASIS     ...  creates a polygonal basis
#  CREATE_POLYNOMIAL_BASIS  ...  creates a polynomial basis
#  CREATE_POWER_BASIS       ...  creates a monomial basis
#  last modified 2007 May 3 by Spencer Graves
#  Previously modified 8 December 2005

if (nargs()==0) {
    type      <- "bspline"
    rangeval  <- c(0,1)
    nbasis    <- 2
    params    <- NULL
    dropind   <- NULL
    quadvals  <- NULL
    values    <- vector("list",0)

    basisobj  <- list(type=type, rangeval=rangeval, nbasis=nbasis, 
                      params=params, dropind=dropind, quadvals=quadvals, 
                      values=values)
    oldClass(basisobj) <- "basisfd"
    basisobj
}

#  if first argument is a basis object, return

if (class(type)=="basisfd"){
    basisobj <- type
    return(basisobj)
}

#  check basistype

type <- use.proper.basis(type)
if (type=="unknown"){
    stop("TYPE unrecognizable.")
}

#  check if QUADVALS is present, and set to default if not

if (missing(quadvals)) quadvals <- NULL
else if(!is.null(quadvals)){
     nquad <- dim(quadvals)[1]
     ncol  <- dim(quadvals)[2]
     if ((nquad == 2) && (ncol > 2)){
         quadvals <- t(quadvals)
         nquad    <- dim(quadvals)[1]
         ncol     <-dim(quadvals)[2]
     }
     if (nquad < 2) stop("Less than two quadrature points are supplied.")
     if (ncol != 2) stop("QUADVALS does not have two columns.")
}

#  check VALUES if present, and set to a single empty cell if not.
if(!missing(values) && !is.null(values)) {
 	 n <- dim(values)[1]
	 k <- dim(values)[2]
    if (n != nquad)   
        stop("Number of rows in VALUES not equal to number of quadrature points.")
    if (k != nbasis)  
        stop("Number of columns in VALUES not equal to number of basis functions.")
}
else values <- NULL

#  check if DROPIND is present, and set to default if not

if(missing(dropind)) dropind<-NULL

#  select the appropriate type and process

if (type=="fourier"){
    paramvec   <- rangeval[2] - rangeval[1]
    period     <- params[1]
    if (period <= 0)  stop("Period must be positive for (a Fourier basis")
    params <- period
    if ((2*floor(nbasis/2)) == nbasis)  nbasis <- nbasis + 1
} else if(type=="bspline"){
    if (!missing(params)){
        nparams  <- length(params)
        if (params[1] <= rangeval[1])       
            stop("Smallest value in BREAKS not within RANGEVAL")
        if (params[nparams] >= rangeval[2]) 
            stop("Largest value in BREAKS not within RANGEVAL")
    }
} else if(type=="expon") {
    if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (exponential basisobj$")
} else if(type=="polyg") {
    if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (polygonal basisobj$")
} else if(type=="power") {
    if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (power basisobj$")
} else if(type=="const") {
    params <- 0
} else if(type=="monom") {
    if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (monomial basisobj$")
} else if(type=="polynom") {
    if (length(params) > 1) 
        stop("More than one parameter for (a polynomial basisobj$")
} else stop("Unrecognizable basis")

#  Save call

obj.call <- match.call()

#  S4 definition
# basisobj <- new("basisfd", call=obj.call, type=type, rangeval=rangeval, nbasis=nbasis, 
#                 params=params, dropind=dropind, quadvals=quadvals, values=values)

#  S3 definition

basisobj <- list(call=obj.call, type=type, rangeval=rangeval, nbasis=nbasis, 
                 params=params, dropind=dropind, quadvals=quadvals, values=values)
oldClass(basisobj) <- "basisfd"

basisobj

}

#  --------------------------------------------------------------------------
#                  print for basisfd class
#  --------------------------------------------------------------------------

print.basisfd <- function(x, ...)
{
  basisobj <- x
  cat("\nBasis object:\n")	
  if (!inherits(basisobj, "basisfd"))
    stop("Argument not a functional data object")
  cat(paste("\n  Type:  ", basisobj$type,"\n"))
  cat(paste("\n  Range: ", basisobj$rangeval[1], 
            " to ",        basisobj$rangeval[2],"\n"))
  if (basisobj$type == "const") return
  cat(paste("\n  Number of basis functions: ", 
            basisobj$nbasis,     "\n"))
  if (basisobj$type == "fourier")
    cat(paste("\n  Period: ",basisobj$params,"\n"))
  if (basisobj$type == "bspline") {
    norder <- basisobj$nbasis - length(basisobj$params)
    cat(paste("\n  Order of spline: ", norder, "\n"))
    if (length(basisobj$params) > 0) {
        print("  Interior knots")
        print(basisobj$params)
    } else {
        print("  There are no interior knots.")
    }
  }
  if (basisobj$type == "polyg") {
    print("  Argument values")
    print(basisobj$params)
  }
  if (basisobj$type == "expon") {
    print("  Rate coefficients")
    print(basisobj$params)
  }
  if (basisobj$type == "power") {
    print("  Exponents")
    print(basisobj$params)
  }
  if (length(basisobj$dropind) > 0) {
    print("  Indices of basis functions to be dropped")
    print(basisobj$dropind)
  }
}

#  --------------------------------------------------------------------------
#                  summary for basisfd class
#  --------------------------------------------------------------------------

summary.basisfd <- function(object, ...)
{
  basisobj <- object
  cat("\nBasis object:\n")	
  if (!inherits(basisobj, "basisfd"))
    stop("Argument not a functional data object")
  cat(paste("\n  Type:  ", basisobj$type,"\n"))
  cat(paste("\n  Range: ", basisobj$rangeval[1], 
            " to ",        basisobj$rangeval[2],"\n"))
  if (basisobj$type == "const") return
  cat(paste("\n  Number of basis functions: ", 
            basisobj$nbasis,     "\n"))
  if (basisobj$type == "fourier")
    cat(paste("\n  Period: ",basisobj$params,"\n"))
  if (length(basisobj$dropind) > 0) {
    print(paste(length(basisobj$dropind),
                "indices of basis functions to be dropped"))
  }
}

#  --------------------------------------------------------------------------
#  equality for basisfd class
#  --------------------------------------------------------------------------

"==.basisfd" <- function(basis1, basis2)
{

# EQ assesses whether two bases are equivalent.

#  Last modified 1 January 2007

type1   <- basis1$type
range1  <- basis1$rangeval
nbasis1 <- basis1$nbasis
pars1   <- basis1$params
drop1   <- basis1$dropind

type2   <- basis2$type
range2  <- basis2$rangeval
nbasis2 <- basis2$nbasis
pars2   <- basis2$params
drop2   <- basis2$dropind

basisequal <- TRUE

#  check types

if (!(type1 == type2)) {
    basisequal <- FALSE
    return(basisequal)
}

#  check ranges

if (range1[1] != range2[1] || range1[2] != range2[2]) {
    basisequal <- FALSE
    return(basisequal)
}

#  check numbers of basis functions

if (nbasis1 != nbasis2) {
    basisequal <- FALSE
    return(basisequal)
}

#  check parameter vectors

if (!(all(pars1 == pars2))) {
    basisequal <- FALSE
    return(basisequal)
}

#  check indices of basis function to drop

if (!(all(drop1 == drop2))) {
    basisequal <- FALSE
    return(basisequal)
}

return(basisequal)

}

#  --------------------------------------------------------------------------
#  pointwise multiplication method for basisfd class
#  --------------------------------------------------------------------------

"*.basisfd" <- function (basisobj1, basisobj2)
{
# TIMES for (two basis objects sets up a basis suitable for (
#  expanding the pointwise product of two functional data
#  objects with these respective bases.  
# In the absence of a true product basis system in this code,
#  the rules followed are inevitably a compromise:
#  (1) if both bases are B-splines, the norder is the sum of the
#      two orders - 1, and the breaks are the union of the
#      two knot sequences, each knot multiplicity being the maximum
#      of the multiplicities of the value in the two break sequences.
#      That is, no knot in the product knot sequence will have a
#      multiplicity greater than the multiplicities of this value
#      in the two knot sequences.  
#      The rationale this rule is that order of differentiability
#      of the product at eachy value will be controlled  by
#      whichever knot sequence has the greater multiplicity.  
#      In the case where one of the splines is order 1, or a step
#      function, the problem is dealt with by replacing the
#      original knot values by multiple values at that location
#      to give a discontinuous derivative.
#  (2) if both bases are Fourier bases, AND the periods are the 
#      the same, the product is a Fourier basis with number of
#      basis functions the sum of the two numbers of basis fns.
#  (3) if only one of the bases is B-spline, the product basis
#      is B-spline with the same knot sequence and order two
#      higher.
#  (4) in all other cases, the product is a B-spline basis with
#      number of basis functions equal to the sum of the two
#      numbers of bases and equally spaced knots.  
 
#  Of course the ranges must also match.
 
#  Last modified 22 March 2007
 
#  check the ranges
 
range1 <- basisobj1$rangeval
range2 <- basisobj2$rangeval
if (range1[1] != range2[1] || range1[2] != range2[2])
    stop("Ranges are not equal.")
 
#  get the types
 
type1 <- basisobj1$type
type2 <- basisobj2$type
 
#  deal with constant bases
 
if (type1 == "const" && type2 == "const") {
    prodbasisobj <- create.constant.basis(range1)
    return(prodbasisobj)
}
 
if (type1 == "const") {
    prodbasisobj <- basisobj2
    return(prodbasisobj)
}
 
if (type2 == "const") {
    prodbasisobj <- basisobj1
    return(prodbasisobj)
}
 
#  get the numbers of basis functions
 
nbasis1 <- basisobj1$nbasis
nbasis2 <- basisobj2$nbasis
 
#  work through the cases
 
if (type1 == "bspline" && type2 == "bspline") {
    #  both are bases B-splines
    #  get orders
    interiorknots1 <- basisobj1$params
    interiorknots2 <- basisobj2$params
    uniqueknots    <- union(interiorknots1, interiorknots2)
    nunique <- length(uniqueknots)
    multunique <- rep(0,nunique)
    for (i in 1:nunique) {
        mult1 <- length(interiorknots1[interiorknots1==uniqueknots[i]])
        mult2 <- length(interiorknots2[interiorknots2==uniqueknots[i]])
        multunique[i] <- max(mult1,mult2)
    }
    allknots <- rep(0,sum(multunique))
    m2 <- 0
    for (i in 1:nunique) {
        m1 <- m2 + 1
        m2 <- m2 + multunique[i]
        allknots[m1:m2] <- uniqueknots[i]
    }
    norder1 <- nbasis1 - length(interiorknots1)
    norder2 <- nbasis2 - length(interiorknots2)
    norder  <- norder1 + norder2 - 1
    allbreaks  <- c(range1[1], allknots, range1[2])
    nbasis <- length(allbreaks) + norder - 2
    prodbasisobj <- 
        create.bspline.basis(range1, nbasis, norder, allbreaks)
    return(prodbasisobj)
}
 
if (type1 == "fourier" && type2 == "fourier") {
    #  both bases Fourier
    #  check whether periods match
    #  if they do not, default to the basis below.
    period1 <- basisobj1$params
    period2 <- basisobj2$params
    nbasis  <- nbasis1 + nbasis2
    if (period1 == period2) {
        prodbasisobj <- create.fourier.basis(range1, nbasis, period1)
        return(prodbasisobj)
    }
}
 
#  default case when all else fails: the product basis is B-spline 
#  When neither basis is a B-spline basis, the order
#  is the sum of numbers of bases, but no more than 8.
#  When one of the bases if B-spline and the other isn"t,
#  the order is the smaller of 8 or the order of the spline
#  plus 2.
 
if (type1 == "bspline" || type2 == "bspline") {
    norder <- 8
    if (type1 == "bspline") {
        interiorknots1 <- basisobj1$params
        norder1        <- nbasis1 - length(interiorknots1)
        norder         <- min(c(norder1+2, norder))
    }
    if (type2 == "bspline") {
        interiorknots2 <- basisobj2$params
        norder2        <- nbasis2 - length(interiorknots2)
        norder         <- min(c(norder2+2, norder))
    }
} else {
    #  neither basis is B-spline
    norder <- min(c(8, nbasis1+nbasis2))
}
#  set up the default B-spline product basis
nbasis <- max(c(nbasis1+nbasis2, norder+1))
prodbasisobj <- create.bspline.basis(range1, nbasis, norder)
return(prodbasisobj)
}
