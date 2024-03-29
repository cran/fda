\name{dirs}
\alias{dirs}
\alias{dirs}
\alias{dirs}
\title{
  Get subdirectories 
}
\description{
  If you want only subfolders and no files, use \code{dirs}.  
  With \code{recursive} = FALSE, \code{\link{dir}} returns both folders
  and files.  With \code{recursive} = TRUE, it returns only files.  
}
\usage{
dirs(path='.', pattern=NULL, exclude=NULL, all.files=FALSE,
     full.names=FALSE, recursive=FALSE, ignore.case=FALSE) 
}
\arguments{
  \item{path, all.files, full.names, recursive, ignore.case}{
    as for \code{\link{dir}}
  }
  \item{pattern, exclude}{
    optional regular expressions of filenames to include or exclude,
    respectively.  
  }
}
\details{
  1.  mainDir <- dir(...)  without recurse 

  2.  Use \code{\link{file.info}} to restrict mainDir to only
  directories.

  3.  If !recursive, return the restricted mainDir.  Else, if
  length(mainDir) > 0, create dirList to hold the results of the
  recursion and call \code{dirs} for each component of mainDir.  Then
  \code{\link{unlist}} and return the result.  
}
\value{
  A character vector of the desired subdirectories.  
}
%\references{}
\author{
  Spencer Graves
} 
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{dir}},
  \code{\link{file.info}}
}
% \examples{
% dirs(path2fdaM)
% dirs(path2fdaM, full.names=TRUE)
% dirs(path2fdaM, recursive=TRUE)
% dirs(path2fdaM, exclude='^@|^private$', recursive=TRUE)

% }
\keyword{IO}% at least one, from doc/KEYWORDS

