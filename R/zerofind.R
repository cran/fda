zerofind <- function(fmat) {
  if (!is.numeric(fmat)) stop("Argument is not numeric.")
  pairNums <- c(min(fmat),max(fmat))
  if (length(pairNums) != 2) stop("Argument is not of length 2.")
  if (min(pairNums) <= 0 && max(pairNums) >= 0) return(TRUE) else return(FALSE)
}
