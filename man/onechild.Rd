\name{onechild}
\alias{onechild}
\title{
  growth in height of one 10-year-old boy 
}
\description{
  Heights of a boy of age approximately 10 collected during one school
  year.  The data were collected "over one school year, with gaps
  corresponding to the school vacations" (AFDA, p. 84)
}
\format{
  A data.frame with two columns:  
    \describe{
	\item{day}{
	  Integers counting the day into data collection with gaps
	  indicating days during which data were not collected.  
	}
	\item{height}{
	  Height of the boy measured on the indicated day.  
	}
    }
}

\source{
  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis, 2nd ed.}, Springer, New York. 

  Tuddenham, R. D., and Snyder, M. M. (1954) "Physical growth of
  California boys and girls from birth to age 18", \emph{University of
  California Publications in Child Development}, 1, 183-364.  
}
%\seealso{\code{\link{eval.penalty}}}
\examples{
with(onechild, plot(day, height, type="b"))

}
% docclass is function
\keyword{datasets}
