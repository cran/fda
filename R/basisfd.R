basisfd <- function(type, rangeval, nbasis, params, dropind=NULL, 
                  quadvals=NULL, values=vector("list",0))
{
#  BASISFD  generator function of "basisfd" class.
#  Arguments:
#  BASISTYPE ... a string indicating the type of basis.  This may be one of
#               "Fourier", "fourier", "Fou", "fou",
#               "Bspline", "bspline", "Bsp", "bsp",
#               "pol", "poly", "polynomial",
#               "mon", "monom", "monomial",
#               "con", "const", "constant"
#               "exp", "exponen", "exponential"
#               "polyg" "polygon", "polygonal"
#  RANGEVAL ... an array of length 2 containing the lower and upper
#               boundaries for the rangeval of argument values
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
#                contains quadrature weights.  For example, for Simpson"s
#                rule for NQUAD = 7, the points are equally spaced and the
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
#  An alternative name for this function is CREATE_BASIS, but PARAMS argument
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

#  last modified 8 December 2005

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
  if(!missing(values) && !is.null(values) && !is.null(quadvals)) {
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
    if (period <= 0)  stop("Period must be positive for a Fourier basis")
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
      stop("No. of parameters not equal to no. of basis fns for exponential basis.")
  } else if(type=="polyg") {
    if (length(params) != nbasis) 
      stop("No. of parameters not equal to no. of basis fns for polygonal basis.")
  } else if(type=="power") {
    if (length(params) != nbasis) 
      stop("No. of parameters not equal to no. of basis fns for power basis.")
  } else if(type=="const") {
    params <- 0
  } else if(type=="monom") {
    if (length(params) != nbasis) 
      stop("No. of parameters not equal to no. of basis fns for monomial basis.")
  } else if(type=="polynom") {
    if (length(params) > 1) 
      stop("More than one parameter for a polynomial basis.")
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

#  "print" method for "basisfd"

print.basisfd <- function(x, ...)
{
	
	cat("\nFunctional data basis x:\n\n")
	
	#  Type
	
	cat(paste(" Type:", x$type), "\n")
	
	#  Range
	
	rangeval <- x$rangeval
	cat(paste(" Range:", rangeval[1], "to", rangeval[2]), "\n")
	
	#  Number of basis functions
	
	cat(paste(" Number of basis functions:", x$nbasis), "\n")
	
   #  display parameters

   cat(" Parameters:\n")
	
   if (x$type == "fourier") {
      cat(paste("  Period: ",x$params),"\n")
   }
   if (x$type == "bspline") {
      norder = x$nbasis - length(x$params)
      cat(paste("  Order of spline: ", norder),     "\n")
      if (length(x$params) > 0) {
        cat("  Interior knots\n")
        print(x$params)
      } else {
        cat("  There are no interior knots.\n")
      }
   }
   if (x$type == "polyg") {
      cat("  Argument values:\n")
      print(x$params)
   }
   if (x$type == "expon") {
      cat("  Rate coefficients:\n")
      print(x$params)
   }
   if (x$type == "power") {
      cat("  Exponents:\n")
      print(x$params)
   }

   #  display indices of basis functions to be dropped

   if (length(x$dropind) > 0) {
      cat("  Indices of basis functions to be dropped:\n")
      print(x$dropind)
   }
	
}


