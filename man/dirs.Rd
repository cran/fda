\name{dirs}
\alias{dirs}
\alias{dirs}
\alias{dirs}
\title{
  Get subdirectories 
}
\description{
  If you want only subfolders and no files,
  use \code{dirs}.  With \code{recursive} =
  FALSE, \code{\link{dir}} returns both
  folders and files.  With \code{recursive} =
  TRUE, it returns only files.  
}
\usage{
dirs(path='.', pattern=NULL, exclude=NULL, 
  all.files=FALSE, full.names=FALSE, 
  recursive=FALSE, ignore.case=FALSE) 
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
\seealso{
  \code{\link{dir}},
  \code{\link{file.info}}
}
\examples{
path2fdaM <- system.file('Matlab/fdaM', package='fda')
dirs(path2fdaM)
dirs(path2fdaM, full.names=TRUE)
dirs(path2fdaM, recursive=TRUE)
dirs(path2fdaM, exclude='^@|^private$', recursive=TRUE)

# Directories to add to Matlab path
# for R.matlab and fda
R.matExt <- system.file('externals', package='R.matlab')
fdaM <- dirs(path2fdaM, exclude='^@|^private$', full.names=TRUE,
              recursive=TRUE)  
add2MatlabPath <- c(R.matExt, path2fdaM, fdaM) 

}
\keyword{IO}% at least one, from doc/KEYWORDS

