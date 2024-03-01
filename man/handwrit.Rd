\name{handwrit}
\docType{data}
\alias{handwrit}
\alias{handwritTime}
\title{
  Cursive handwriting samples
}
\description{
  20 cursive samples of 1401 (x, y,) coordinates for writing "fda"
}
\usage{
handwrit
handwritTime
}
\format{
  \describe{
    \item{handwrit}{
      An array of dimensions (1401, 20, 2) giving 1401 pairs of (x, y)
      coordinates for each of 20 replicates of cursively writing "fda"
    }
    \item{handwritTime}{
      seq(0, 2300, length=1401) = sampling times
    }
  }
}
\details{
  These data are the X-Y coordinates of 20 replications of writing
  the script "fda".  The subject was Jim Ramsay.  Each replication
  is represented by 1401 coordinate values.  The scripts have been
  extensively pre-processed.  They have been adjusted to a common
  length that corresponds to 2.3 seconds or 2300 milliseconds, and
  they have already been registered so that important features in
  each script are aligned.

  This analysis is designed to illustrate techniques for working
  with functional data having rather high frequency variation and
  represented by thousands of data points per record.  Comments
  along the way explain the choices of analysis that were made.

  The final result of the analysis is a third order linear
  differential equation for each coordinate forced by a
  constant and by time.  The equations are able to reconstruct
  the scripts to a fairly high level of accuracy, and are also
  able to accommodate a substantial amount of the variation in
  the observed scripts across replications.  by contrast, a
  second order equation was found to be completely inadequate.

  An interesting surprise in the results is the role placed by
  a 120 millisecond cycle such that sharp features such as cusps
  correspond closely to this period.  This 110-120 msec cycle
  seems is usually seen in human movement data involving rapid
  movements, such as speech, juggling and so on.

  These 20 records have already been normalized to a common time
  interval of 2300 milliseconds and have beeen also registered so that
  prominent features occur at the same times across replications.  Time
  will be measured in (approximate) milliseconds and space in meters.
  The data will require a small amount of smoothing, since an error of
  0.5 mm is characteristic of the OPTOTRAK 3D measurement system used to
  collect the data.

  Milliseconds were chosen as a time scale in order to make the ratio of
  the time unit to the inter-knot interval not too far from one.
  Otherwise, smoothing parameter values may be extremely small or
  extremely large.

  The basis functions will be B-splines, with a spline placed at each
  knot.  One may question whether so many basis functions are required,
  but this decision is found to be essential for stable derivative
  estimation up to the third order at and near the boundaries.

  Order 7 was used to get a smooth third derivative, which requires
  penalizing the size of the 5th derivative, which in turn requires an
  order of at least 7.  This implies norder + no. of interior knots =
  1399 + 7 = 1406 basis functions.

  The smoothing parameter value 1e8 was chosen to obtain a fitting error
  of about 0.5 mm, the known error level in the OPTOTRACK equipment.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\examples{
oldpar <- par(no.readonly=TRUE)
plot(handwrit[, 1, 1], handwrit[, 1, 2], type="l")
par(oldpar)
}
\keyword{datasets}

