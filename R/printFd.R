printFd <- function(fd, ...)
{
  #  Prints a functional data object FD
  #  The remaining optional arguments are the same as those available
  #     in the regular "print" function.


 #  Last modified 23 March 2003

  if (inherits(fd, "fd")) {
    cat("Dimensions of data:\n")
    print(fd$fdnames)
  } else {
    stop("First argument is not a functional data object.")
  }
  cat("Coefficient Matrix:\n")
  coef <- getcoef(fd)
  print(coef)
  fbdo   <- getbasis(fd)
  type   <- fbdo$type
  params <- fbdo$params
  cat("\nBasis:\n")
  cat(paste("  Type:", type,"\n"))
  cat(paste("  Range:",fbdo$rangeval[1],"to",fbdo$rangeval[2],"\n"))
  cat(paste("  Number of basis functions:",fbdo$nbasis,"\n"))
  if (type == "fourier") cat("  Period: ")
  if (type == "bspline") cat("  Interior knots         \n")
  if (type == "expon")   cat("  Rate coefficients      \n")
  if (type == "poly")    cat("  Polynomial coefficients\n")
  if (type == "polyg")   cat("  Argument values        \n")
  cat(format(params))
  cat("\n")
}
