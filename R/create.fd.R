create.fd <- function (coef, basisfd, fdnames=defaultnames)
{

  #  This function creates a functional data object.
  #    A functional data object consists of a basis for expanding a functional
  #    observation and a set of coefficients defining this expansion.
  #    The basis is contained in a "basis.fd" object; that is, a realization
  #    of the "basis.fd" class.

  #  Arguments
  #  COEF ... An array containing coefficient values for the expansion of each
  #             set of function values in terms of a set of basis functions.
  #           If COEF is a three-way array, then the first dimension
  #             corresponds to basis functions, the second to replications,
  #             and the third to variables.
  #           If COEF is a matrix, it is assumed that there is only
  #             one variable per replication, and then
  #                 rows    correspond to basis functions
  #                 columns correspond to replications
  #           If COEF is a vector, it is assumed that there is only one
  #             replication and one variable.
  #  BASISFD ... a functional data basis object
  #  FDNAMES  ... The analogue of the dimnames attribute of an array, this is
  #               a list of length 3 with members containing:
  #               1. a character vector of names for the argument values
  #               2. a character vector of names for the replications or cases
  #               3. a character vector of names for the functions
  #               Each of these vectors can have a name referring to the modality
  #                 of the data.  An example would be "time", "reps", "values"

  #  Returns:
  #  FD ... a functional data object

  #  Note:  Earlier versions also had members named "df" and "gcv".  These
  #         have been removed.

  #  last modified 29 April 2003

  #  check COEF and get its dimensions

  if(!is.numeric(coef)) stop("coef must be numerical vector or matrix")
  else if (is.vector(coef)) {
            coef  <- as.matrix(coef)
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
  else if (is.matrix(coef)) {
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
  else if (is.array(coef)) {
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
  else stop("argument coef is not correct")

  if (ndim > 3) stop(
      "First argument not of dimension 1, 2 or 3")

  #  check BASISFD

  if (!(inherits(basisfd, "basis.fd"))) stop(
    "Argument BASISFD must be of basis.fd class")

  if (dim(coef)[1] != basisfd$nbasis)
    stop("Number of coefficients does not match number of basis functions.")

  #  setup number of replicates and number of variables

  if (ndim > 1) nrep <- coefd[2] else nrep <- 1
  if (ndim > 2) nvar <- coefd[3] else nvar <- 1

  #  set up default fdnames

  if (ndim == 1) {
    defaultnames <- list("time", "reps", "values")
  }
  if (ndim == 2) {
    defaultnames <- list("time",
                         paste("reps",as.character(1:nrep)),
                         "values")
  }
  if (ndim == 3) {
    defaultnames <- list("time",
                         paste("reps",as.character(1:nrep)),
                         paste("values",as.character(1:nvar)) )
  }
  names(defaultnames) <- c("args", "reps", "funs")

  fd        <- list( coef, basisfd, fdnames )
  names(fd) <- c("coefs", "basis", "fdnames")
  class(fd) <- "fd"

  return(fd)
}
