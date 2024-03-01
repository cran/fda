\name{plotscores}
\alias{plotscores}
\title{
  Plot Principal Component Scores
}
\description{
The coefficients multiplying the harmonics or principal component functions
are plotted as points.
}
\usage{
plotscores(pcafd, scores=c(1, 2), xlab=NULL, ylab=NULL,
           loc=1, matplt2=FALSE, ...)
}
\arguments{
\item{pcafd}{
an object of the "pca.fd" class that is output by function
\code{pca.fd}.
}
\item{scores}{
the indices of the harmonics for which coefficients are
plotted.
}
\item{xlab}{
a label for the horizontal axis.
}
\item{ylab}{
a label for the vertical axis.
}
\item{loc}{
an integer:
if loc  >0, you can then click on the plot in loc places and you'll get
plots of the functions with these values of the principal component
coefficients.
}
\item{matplt2}{
a logical value:
if \code{TRUE}, the curves are plotted on the same plot;
otherwise, they are plotted separately.
}
\item{\dots }{
additional plotting arguments used in function \code{plot}.
}
}
\section{Side Effects}{
a plot of scores
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
  \code{\link{pca.fd}}
}
\keyword{smooth}
