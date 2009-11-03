\name{fdaMatlabPath}
\alias{fdaMatlabPath}
\title{
  Add 'fdaM' to the Matlab path 
}
\description{
  Write a sequence of Matlab commands to \code{fdaMatlabPath.m} in the
  working directory containing commands to add \code{fdaM} to the path
  for Matlab.  
}
\usage{
fdaMatlabPath(R.matlab) 
}
\arguments{
  \item{R.matlab}{
    logical:  If TRUE, include '~R/library/R.matlab/externals' in the
    path.  If(missing(R.matlab)) include '~R/library/R.matlab/externals'
    only if R.matlab is installed.   
  }
}
\details{
  \itemize{
    \item{USAGE}{
      If your Matlab installation does NOT have a \code{startup.m} file,
      it might be wise to copy \code{fdaMatlabPath.m} into a directory
      where Matlab would look for \code{startup.m}, then rename it to
      \code{startup.m}.

      If you have a \code{startup.m}, you could add the contents of
      \code{fdaMatlabPath.m} to \code{startup.m}.

      Alternatively, you can copy \code{fdaMatlabPath.m} into the
      directory containing \code{startup.m} and add the following to the
      end of \code{startup.m}:

      \tabular{lll}{
	\tab if exist('fdaMatlabPath') \tab \cr
	\tab \tab fdaMatlabPath ; \cr
	\tab end \tab 
      }
    }
    \item{ALGORITHM}{
      1.  path2fdaM = path to the \code{Matlab/fdaM} subdiretory of the
      \code{fda} installation directory. 

      2.  Find all subdirectories of path2fdaM except those beginning in
      '@' or including 'private'. 

      3.  if(requires(R.matlab)) add the path to \code{MatlabServer.m}
      to \code{dirs2add} 

      4.  d2a <- paste("addpath('", dirs2add, "');", sep='')

      5.  writeLines(d2a, 'fdaMatlabPath.m')

      6.  if(exists(startupFile)) append \code{d2a} to it
    }
  }
}
\value{
  A character vector of Matlab \code{addpath} commands is returned
  invisibly.  
}
\references{
  Matlab documentation for \code{addpath} and \code{startup.m}.  
}
\author{
  Spencer Graves with help from Jerome Besnard
} 
\seealso{
  \code{\link[R.matlab]{Matlab}},
  \code{\link{dirs}}
}
\examples{
# Modify the Matlab startup.m only when you really want to,
# typically once per installation ... certaintly not
# every time we test this package.
fdaMatlabPath()
}
\keyword{programming}% at least one, from doc/KEYWORDS

