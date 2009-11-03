\name{file.copy2}
\alias{file.copy2}
\title{
  Copy a file with a default 'to' name 
}
\description{
  Copy a file appending a number to make the \code{to} name unique, with
  default \code{to} = \code{from}.    
}
\usage{
file.copy2(from, to) 
}
\arguments{
  \item{from}{ character:  name of a file to be copied }
  \item{to}{   
    character:  name of copy.  Default = \code{from} with an integer
    appended to the name.  
  }
}
\details{
  1.  length(from) != 1:  Error:  Only one file can be copied.  
  
  2.  file.exists(from)?  If no, If no, return FALSE.

  3.  if(missing(to))to <- from;  else if(length(to)!=1) error.

  4.  file.exists(to)?  If yes, Dir <- dir(dirname(to)), find all
  \code{Dir} starting with \code{to}, and find the smallest integer to
  append to make a unique \code{to} name.  

  5.  file.copy(from, to)

  6.  Return TRUE.  
}
\value{
  logical:  TRUE (with a name = name of the file created);  FALSE if no
  file created.  
}
%\references{}
\author{
  Spencer Graves 
} 
\seealso{
  \code{\link{file.copy}},
%  \code{\link{fdaMatlabPath}}
}
\examples{
\dontrun{
file.copy2('startup.m')
# Used by 'fdaMatlabPath' so an existing 'startup.m' is not destroyed
}
}
\keyword{IO}% at least one, from doc/KEYWORDS

