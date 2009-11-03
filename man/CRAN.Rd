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
CRAN(CRAN_pattern, n_R_CHECK4CRAN)
}
\arguments{
  \item{CRAN_pattern}{
    a regular expressions to apply to the names of \code{Sys.getenv()}
    to identify possible CRAN parameters.  Defaults to
    \code{Sys.getenv('_CRAN_pattern_')} if available and '^_R_' if not.
  }
  \item{n_R_CHECK4CRAN}{
    Assume this is CRAN if at least n_R_CHECK4CRAN elements of
    \code{Sys.getenv()} have names matching x.  Defaults to
    \code{Sys.getenv('_n_R_CHECK4CRAN_')} if available and 5 if not.
  }
}
\details{
  The "Writing R Extensions" manual says that "R CMD check" can be
  customized "by setting environment variables _R_CHECK_*_:, as
  described in" the Tools section of the "R Internals" manual.

  'R CMD check' was tested with R 3.0.1 under Fedora 18 Linux and with
  Rtools 3.0 from April 16, 2013 under Windows 7.  With the '--as-cran'
  option, 7 matches were found;  without it, only 3 were found.  These
  numbers were unaffected by the presence or absence of the '--timings'
  parameter.  On this basis, the default value of n_R_CHECK4CRAN was set
  at 5.

  1.  x. <- Sys.getenv()

  2.  Fix \code{CRAN_pattern} and \code{n_R_CHECK4CRAN} if missing.

  3.  Let i be the indices of x. whose names match all the patterns in
  the vector x.

  4.  Assume this is CRAN if length(i) >= n_R_CHECK4CRAN

}
\value{
  a logical scalar with attributes 'Sys.getenv' containing the results
  of \code{Sys.getenv()} and 'matches' contining \code{i} per step 3
  above.
}

\seealso{
  \code{\link{Sys.getenv}}
}
\examples{
% Look in fda-ex.Rout
% for the results from R CMD check
% Modified defaults were tested with environment variables
% _CRAN_pattern_ = 'A' and _n_R_CHECK4CRAN_ = '10'
% 'R CMD check' found 26 matches and CRAN() returned TRUE
cran <- CRAN()
str(cran)
gete <- attr(cran, 'Sys.getenv')
(ngete <- names(gete))

iget <- grep('^_', names(gete))
gete[iget]

\dontrun{
if(CRAN()){
  stop('CRAN')
} else {
  stop('NOT CRAN')
}
}
}
\keyword{environment}

