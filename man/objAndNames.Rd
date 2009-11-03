\name{objAndNames}
\alias{objAndNames}
\title{
  Add names to an object 
}
\description{
  Add names to an object from 'preferred' if available and 'default' if
  not.  
}
\usage{
objAndNames(object, preferred, default)
}
\arguments{
  \item{object}{
    an object of some type to which names must be added.  If
    length(dim(object))>0 add 'dimnames', else add 'names'.  
  }
  \item{preferred}{
    A list to check first for names to add to 'object'.  
  }
  \item{default}{
    A list to check for names to add to 'object' if appropriate names
    are not found in 'preferred'.  
  }
}
\details{
  1.  If length(dim(object))<1, names(object) are taken from 'preferred'
  if they are not NULL and have the correct length, else try 'default'.

  2.  Else for(lvl in 1:length(dim(object))) take dimnames[[lvl]] from
  'preferred[[i]]' if they are not NULL and have the correct length,
  else try 'default[[lvl]].  
}
\value{
  An object of the same class and structure as 'object' but with either
  names or dimnames added or changed.  
}
\author{ Spencer Graves  }

\seealso{
  \code{\link{bifd}}
}
\examples{
# The following should NOT check 'anything' here
tst1 <- objAndNames(1:2, list(letters[1:2], LETTERS[1:2]), anything)
all.equal(tst1, c(a=1, b=2))

# The following should return 'object unchanged
tst2 <- objAndNames(1:2, NULL, list(letters))
all.equal(tst2, 1:2)


tst3 <- objAndNames(1:2, list("a", 2), list(letters[1:2]))
all.equal(tst3, c(a=1, b=2) )

# The following checks a matrix / array
tst4 <- array(1:6, dim=c(2,3))
tst4a <- tst4
dimnames(tst4a) <- list(letters[1:2], LETTERS[1:3])
tst4b <- objAndNames(tst4, 
       list(letters[1:2], LETTERS[1:3]), anything)
all.equal(tst4b, tst4a)

tst4c <- objAndNames(tst4, NULL,        
       list(letters[1:2], LETTERS[1:3]) )
all.equal(tst4c, tst4a)

}
% docclass is function
\keyword{attribute}
