summaryFd <- function(fd,...)
{
  #  Summarizes a functional data object FD or a basis object
  #  The remaining optional arguments are the same as those available
  #     in the regular "summary" function.

  #  Last modified 6 Feb 2001

  if (inherits(fd, "fd")) {
    cat("Dimensions of data:\n")
    print(fd$fdnames)
  } else {
    stop("First argument is not a functional data object.")
  }

  fbdo <- getbasis(fd)
  cat("\nBasis:\n")
  cat(paste("  Type:", fbdo$type,"\n"))
  cat(paste("  Range:",fbdo$rangeval[1],"to",fbdo$rangeval[2],"\n"))
  cat(paste("  Number of basis functions:",  fbdo$nbasis,     "\n"))
  type <- getbasistype(fbdo)
  if  (type == "fourier") {
    cat(paste("  Period:",fbdo$params,"\n"))
  }
  if (type == "bspline") {
    cat("  Interior knots\n")
    print(fbdo$params)
  }
  if (type == "poly") {
    cat("  Polynomial coefficients\n")
    print(fbdo$params)
  }
  if (type == "polyg") {
    cat("  Argument values\n")
    print(fbdo$params)
  }
  if (type == "expon") {
    cat("  Rate coefficients\n")
    print(fd[[2]]$params)
  }
}
