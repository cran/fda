#  "summary" method for "fd"
summary.fd <- function(object,...)
{
	
	cat("Functional data object:\n\n")
	
	cat(" Dimensions of the data:\n")
	cat(paste("   ",names(object$fdnames),"\n"))
	
	print.basisfd(object$basis)
	
	cat("\nCoefficient matrix:\n\n")
	
	object$coefs
	
}