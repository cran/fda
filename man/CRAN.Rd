\name{CRAN}
\alias{CRAN}
\title{
  Test if running as CRAN
}
\description{
  This function allows package developers to run tests themselves that
  should not run on CRAN or with "R CMD check --as-cran" because of
  compute time constraints with CRAN tests.
}
\usage{
CRAN(x='_R_CHECK_TIMINGS_')
}
\arguments{
  \item{x}{
    name of an environmental variable to check.
  }
}
\details{
  x. <- Sys.getenv(x)
  xl <- as.logical(x.)

  return((!is.na(xl)) || xl)

  This returns FALSE either if the environmental variable is not found
  or it is not TRUE.

  Duncan Murdoch noted that the following variables are set by
  "--as-cran":  _R_CHECK_TIMINGS_, _R_CHECK_INSTALL_DEPENDS_,
  _R_CHECK_NO_RECOMMENDED_, and _R_SHLIB_BUILD_OBJECTS_SYMBOL_TABLES_.
  For this purpose, he suggested we use the first.
}
\value{
  a logical scalar
}

%\seealso{}
\examples{
\dontrun{
cran <- CRAN()
str(cran)
print(cran)
stop()

if(CRAN()){
  stop('CRAN')
} else {
  stop('NOT CRAN')
}
}
}
\keyword{environment}

