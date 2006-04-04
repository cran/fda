\name{fd}
\alias{fd}
\title{
  Define a Functional Data Object
}
\description{
This is the constructor function for objects of the \code{fd} class.
Each function that sets up an object of this class must call this function.
This includes functions \code{data2fd},
\code{smooth.basis}, \code{density.fd},  and
so forth that estimate functional data objects that smooth or otherwise
represent data.  Ordinarily, user
of the functional data analysis software will not need to call this function
directly, but these notes are valuable to understanding what the "slots"
or "members" of the \code{fd} class are.
}
\usage{
fd(coef=matrix(0, 2, 1), basisobj=basisfd(),
   fdnames=defaultnames)
}
\arguments{
\item{coef}{
a vector, matrix, or three-dimensional array of coefficients.
The first dimension of
matrices and arrays, as well as that of a vector, is corresponds to
basis functions.  The second dimension corresponds to the number of
functions or curves, or to replicates.  The third dimension if present
corresponds to variables for multivariate functional data objects.
A functional data object is referred to as "univariate" if this
argument is a vector or a matrix, and "multivariate" if it is a
three-dimensional array.
}
\item{basisobj}{
a functional basis object defining the basis
}
\item{fdnames}{
A list of length 3, each member being a string vector containing labels for the
levels of the corresponding dimension of the discrete data.  The first dimension
is for argument values, and is given the default name "time", the second is
for replications, and is given the default name "reps", and the
third is for
functions, and is given the default name "values".
}
}
\value{
a functional data object
}
\details{


 To check that an object is of this class, use function
\code{is.fd}.

 Normally only developers of new functional data analysis
functions will actually need
to use this function, so no examples are provided.
}
\seealso{
Most of the functional data analysis functions either have
functional data objects as arguments, or generate functional data
objects, or both.
}

% docclass is function
\keyword{smooth}
