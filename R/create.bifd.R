create.bifd <- function (coef, sbasisfd, tbasisfd,
                         bifdnames = list(NULL, repnames, NULL ))
{

  #  This function creates a bi-functional data object.  A bi-functional datum
  #  object consists of two bases for expanding a bivariate function and
  #  a set of coefficients defining this expansion.  Each basis is contained
  #  in a "basis.fd" object.  That is, a realization of the "basis.fd" class.

  #  Arguments
  #  COEF     ... a two-, three-, or four-dimensional array containing
  #               coefficient values for the expansion of each set of bivariate
  #               function values in terms of a set of basis function values
  #               If COEF is a two-way, it is assumed that there is only
  #                 one variable and only one replication, and then
  #                 the first and second dimensions correspond to
  #                 the basis functions for the first and second argument,
  #                 respectively.
  #               If COEF is a three-way, it is assumed that there is only
  #                 one variable per replication, and then
  #                 the first and second dimensions correspond to
  #                 the basis functions for the first and second argument,
  #                 respectively, and the third dimension corresponds to
  #                 replications.
  #               If COEF is a four-way array, then the fourth dimension
  #                 corresponds to variables
  #  SBASISFD ... a functional data basis object for the first  argument s
  #  TBASISFD ... a functional data basis object for the second argument t
  #  BIFDNAMES ... A list of length 3 with members containing
  #               1. a single name for the argument domain, such as "Time"
  #               2. a vector of names for the replications or cases
  #               3. a name for the function, or a vector of names if there
  #                  are multiple functions.
  #  Returns
  #  BIFDOBS  ... a functional datum object


 #  Last modified 6 Feb 2001

  if (length(dim(coef)) == 2) {
    repnames <- NULL
  } else {
    repnames <- dimnames(coef)[3]
  }

  if (is.vector(coef) || length(dim(coef)) > 4) stop(
      " First argument not of dimension 2, 3 or 4")

  if (!(inherits(sbasisfd, "basis.fd"))) stop(
    "Argument SBASISFD must be of basis.fd class")
  if (!(inherits(tbasisfd, "basis.fd"))) stop(
    "Argument TBASISFD must be of basis.fd class")
  if (dim(coef)[1] != sbasisfd$nbasis) stop(paste(
         "First dimension does not match number of basis functions",
         "for first argument."))
  if (dim(coef)[2] != tbasisfd$nbasis) stop(paste(
         "Second dimension does not match number of basis functions",
         "for second argument."))

  bifd        <- list( coef, sbasisfd, tbasisfd, bifdnames)
  names(bifd) <- c("coefs", "sbasis", "tbasis", "bifdnames")
  setOldClass("bifd")
  oldClass(bifd) <- "bifd"

  return(bifd)
}
