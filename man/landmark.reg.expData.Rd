\name{landmark.reg.expData}
\docType{data}
\alias{landmark.reg.expData}
\alias{sampleData}

\title{
  Experiment data for landmark registration and alignment
}
\description{
  A data.frame object obtaining flow cytometry data with 5 samples
  subject to landmark alignment. 
}

\format{
  \describe{
    \item{data}{Flow cytometry fluorescent intensity data.}
    \item{which}{Indicating the sample name of the data.}
  }  
}

\usage{
sampleData
}

\seealso{\code{\link{landmarkreg}}}
\examples{
data(landmark.reg.expData)
head(sampleData)
}

\author{C. J. Wong \email{cwon2@fhcrc.org}}

% docclass is function
\keyword{datasets}
