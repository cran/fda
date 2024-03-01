\name{pcaPACE}
\alias{pcaPACE}
\title{
  Estimate the functional principal components 
}
\description{
  Carries out a functional PCA with regularization from the estimate of the covariance surface
}
\usage{
  pcaPACE(covestimate, nharm, harmfdPar, cross)
}
\arguments{
  \item{covestimate}{
    a list with the two named entries "cov.estimate" and "meanfd"
  }
  \item{nharm}{
    the number of harmonics or principal components to compute.
  }
  \item{harmfdPar}{
    a functional parameter object that defines the
    harmonic or principal component functions to be estimated.
  }
  \item{cross}{
    a logical value:
    if TRUE, take into account the cross covariance for estimating the eigen functions.
  }
}
\value{
  an object of class "pca.fd" with these named entries:
    
    \item{harmonics}{
      a functional data object for the harmonics or eigenfunctions
    }
  \item{values}{
    the complete set of eigenvalues
  }
  \item{scores}{
    NULL. Use "scoresPACE" for estimating the pca scores
  }
  \item{varprop}{
    a vector giving the proportion of variance explained
    by each eigenfunction
  }
  \item{meanfd}{
    a functional data object giving the mean function
  }
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.

  Yao, F., Mueller, H.G., Wang, J.L. (2005), 
  \emph{Functional data analysis for sparse longitudinal data},
  J. American Statistical Association, 100, 577-590.
}
