wtcheck = function(n, wtvec) {

if (is.null(wtvec)) wtvec = rep(1,n)
if (!is.vector(wtvec))  stop("'wtvec' is not a vector.")
if (length(wtvec) != n) stop("'wtvec' of wrong length")
if (min(wtvec) <= 0)    stop("All values of 'wtvec' must be positive.")

return(wtvec)

}
