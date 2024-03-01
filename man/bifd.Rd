\name{bifd}
\alias{bifd}
\title{
  Create a bivariate functional data object 
}
\description{
  This function creates a bivariate functional data object, which
  consists of two bases for expanding a functional data object of two 
  variables, s and t, and a set of coefficients defining this expansion.
  The bases are contained in "basisfd" objects.
}
\usage{
bifd (coef=matrix(0,2,1), sbasisobj=create.bspline.basis(),
      tbasisobj=create.bspline.basis(), fdnames=defaultnames)
}
\arguments{
  \item{coef}{
    a two-, three-, or four-dimensional array containing
    coefficient values for the expansion of each set of bivariate
    function values=terms of a set of basis function values

    If 'coef' is two dimensional, this implies that there is only
    one variable and only one replication.  In that case, 
    the first and second dimensions correspond to
    the basis functions for the first and second argument,
    respectively.
    
    If 'coef' is three dimensional, this implies that there are multiple
    replicates on only one variable.  In that case, 
    the first and second dimensions correspond to
    the basis functions for the first and second argument,
    respectively, and the third dimension corresponds to
    replications.
    
    If 'coef' has four dimensions, the fourth dimension
    corresponds to variables.
  }
  \item{sbasisobj}{
    a functional data basis object
    for the first argument s of the bivariate function.  
  }
  \item{tbasisobj}{
    a functional data basis object
    for the second argument t of the bivariate function.  
  }
  \item{fdnames}{
    A list of length 4 containing dimnames for 'coefs' if it is a
    4-dimensional array.  If it is only 2- or 3-dimensional, the later
    components of fdnames are not applied to 'coefs'.  In any event, the
    components of fdnames describe the following:

    (1) The row of 'coefs' corresponding to the bases in sbasisobj.
    Defaults to sbasisobj[["names"]] if non-null and of the proper
    length, or to existing dimnames(coefs)[[1]] if non-null and of the
    proper length, and to 's1', 's2', ...,
    otherwise.  

    (2) The columns of 'coefs' corresponding to the bases in tbasisobj.
    Defaults to tbasisobj[["names"]] if non-null and of the proper
    length, or to existing dimnames(coefs)[[2]] if non-null and of the
    proper length, and to 't1', 't2', ...,  
    otherwise.  

    (3) The replicates.  Defaults to dimnames(coefs)[[3]] if non-null
    and of the proper length, and to 'rep1', ..., otherwise.

    (4) Variable names.  Defaults to dimnames(coefs)[[4]] if non-null
    and of the proper length, and to 'var1', ..., otherwise.  
  }
}

\value{
  A bivariate functional data object = a list of class 'bifd' 
  with the following components:

  \item{coefs}{
    the input 'coefs' possible with dimnames from dfnames if provided or
    from sbasisobj$names and tbasisobsj$names
  }
  \item{sbasisobj}{
    a functional data basis object
    for the first argument s of the bivariate function.  
  }
  \item{tbasisobj}{
    a functional data basis object
    for the second argument t of the bivariate function.  
  }
  \item{bifdnames}{
    a list of length 4 giving names for the dimensions of coefs, with
    one or two unused lists of names if length(dim(coefs)) is only two
    or one, respectively.  
  }
}
\author{ Spencer Graves  }

\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
  \emph{Functional data analysis with R and Matlab}, Springer, New
  York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{basisfd}}
  \code{\link{objAndNames}}
}
\examples{
Bspl2 <- create.bspline.basis(nbasis=2, norder=1)
Bspl3 <- create.bspline.basis(nbasis=3, norder=2)

(bBspl2.3 <- bifd(array(1:6, dim=2:3), Bspl2, Bspl3))
str(bBspl2.3)

}
% docclass is function
\keyword{attribute}
