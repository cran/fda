create.basis.fd <- function (type, rangeval, nbasis, params = paramvec)
{

  #  This function creates a functional data basis.
  #  Arguments
  #  TYPE     ... a string indicating the type of basis.  This may be one of
  #               "Fourier", "fourier", "Fou", "fou",
  #               "Bspline", "bspline", "Bsp", "bsp",
  #               "Polynomial","polynomial","Polyn", "polyn",
  #               "Exponential","exponential","Exp", "exp",
  #               "Polygonal", "polygonal","polyg" "polyg",
  #               "Constant","constant","Con", "con",
  #               "power","power","Pow", "pow",
  #  RANGEVAL ... an array of length 2 containing the lower and upper
  #               boundaries for the rangeval of argument values
  #  NBASIS   ... the number of basis functions
  #  PARAMS   ... If the basis is "fourier", this is a single number indicating
  #                 the period.  That is, the basis functions are periodic on
  #                 the interval [0,PARAMS] or any translation of it.
  #               If the basis is "bspline", the values are interior points at
  #                 which the piecewise polynomials join.
  #                 Note that the number of basis functions NBASIS is equal
  #                 to the order of the Bspline functions plus the number of
  #                 interior knots, that is the length of PARAMS.
  #                 This means that NBASIS must be at least 1 larger than the
  #                 length of PARAMS.
  #               If the basis is "polynomial", this is a single number "ctr"
  #                 indicating the zero value for the polynomials, which are
  #                 of the form (x - ctr)^m, m=0,...,NBASIS-1.
  #               If the basis is "exponential", this is a vector of rate
  #                 constants, and the basis functions are of the form
  #                 exp(rate*x)
  #               If the basis is "polygonal" or "constant", value(s) in PARAMS
  #                 are not used.
  #  Returns
  #  BASIS.FD  ... a functional data basis object
  #  An alternative name for this function is CREATE.BASIS, but PARAMS argument
  #     must be supplied.
  #  Specific types of bases may be set up more conveniently using functions
  #  CREATE.BSPLINE.BASIS  ...  creates a b-spline basis
  #  CREATE.FOURIER.BASIS  ...  creates a fourier basis
  #  CREATE.POLYGONal.BASIS  ...  creates a polygonal basis
  #  CREATE.POLYnomial.BASIS  ...  creates a polynomial basis
  #  CREATE.exponential.BASIS  ...  creates a exponential basis
  #  CREATE.power.BASIS  ...  creates a power basis
  #  CREATE.const.BASIS  ...  creates a const basis

  #  Last modified 25 March 2003

  #  Check TYPE

  type <- use.proper.basis(type)
  if (type == "unknown") stop ("TYPE unrecognizable.")

  #  check RANGE

  if (!rangechk(rangeval)) stop('Argument RANGEVAL is not correct.')

  #  check NBASIS

  if (!is.numeric(nbasis) || nbasis<=0)
    stop('nbasis must be a positive integer')
  nbasis <- ceiling(nbasis)

  #  Set up basis depending on type

   if (type == "fourier") {
     paramvec   <- rangeval[2] - rangeval[1]   # default value of period
     period     <- params[1]
     basisfd <- create.fourier.basis(rangeval,nbasis,period)
   }

   if (type == "bspline") {
     basisfd <- create.bspline.basis(rangeval,nbasis,4,params)
   }

   if (type == "polynomial") {
     basisfd <- create.polynomial.basis(rangeval,nbasis,params[1])
   }

   if (type == "exponential") {
     basisfd <- create.exponential.basis(rangeval,params)
   }

   if (type == "polygonal") {
     basisfd <- create.polygonal.basis(params)
   }

   if (type == "constant") {
     basisfd <- create.constant.basis(rangeval)
   }

   if (type == "power") {
     basisfd <- create.power.basis(rangeval,params)
   }

   return(basisfd)
}
