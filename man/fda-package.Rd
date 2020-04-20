\name{fda-package}
\alias{fda-package}
\alias{fda}
\docType{package}
\title{Functional Data Analysis in R}
\description{
  Functions and data sets companion to Ramsay, J. O.; Hooker, Giles; and
  Graves, Spencer (2010) Functional Data Analysis with R and Matlab,
  plus Ramsay, J. O., and Silverman, B. W. (2006) Functional Data
  Analysis, 2nd ed. and (2002) Applied Functional Data Analysis
  (Springer).  This includes finite bases approximations (such as
  splines and Fourier series) to functions fit to data smoothing on the
  integral of the squared deviations from an arbitrary differential operator.
}
\details{
  \tabular{ll}{
    Package: \tab fda\cr
    Type: \tab Package\cr
    Version: \tab 2.2.6\cr
    Date: \tab 2011-02-03\cr
    License: \tab GPL-2\cr
    LazyLoad: \tab yes\cr
  }
}
\author{
  J. O. Ramsay,  <ramsay@psych.mcgill.ca>, Hadley Wickham
  <h.wickham@gmail.com>, Spencer Graves <spencer.graves@prodsyse.com>,
  Giles Hooker <gjh27@cornell.edu>

  Maintainer:  J. O. Ramsay <ramsay@psych.mcgill.ca>
}
\references{
  Ramsay, J. O.; Hooker, Giles; and Graves, Spencer (2010)
  \emph{Functional Data Analysis with R and Matlab}, Springer, New
  York.

  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\examples{
##
## As noted in the Preface to Ramsay, Hooker and Graves (p. v),
## the fda package includes scripts to reproduce all but one of the
## figures in the book.
##
## These figures can be found and run as follows:
##
\dontrun{
scriptsDir <- system.file('scripts', package='fda')
Rscripts <- dir(scriptsDir, full.names=TRUE, pattern='R$')
fdarm <- grep('fdarm', Rscripts, value=TRUE)
chapters <- length(fdarm)
# NOTE:  If R fails in any of these scripts,
# this for loop will not end normally,
# and the abnormal termination will be displayed:
for(ch in 1:chapters){
  cat('Running', fdarm[ch], '\n')
  invisible(source(fdarm[ch]))
}
}
}
\keyword{smooth}
