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
