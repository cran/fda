\name{NEWS}
\alias{NEWS}
\title{Recent changes to the fda package}
\description{
  \itemize{
    \item{Changes in version fda_6.0.3 2022-05-02:}{
      \itemize{
        \item{Landmark registration:}{Landmark registration using function 
          \code{landmarkreg} can no longer be done by using
          function \code{smooth.basis} instead of function \code{smooth.morph}.  The 
          warping function must be strictly monotonic, and we have found that using 
          \code{smooth.basis} too often violates this monotonicity constraint.  Function 
          \code{smooth.morph} ensures monotonicity and in most applications takes negligible
          computer time to do so.
        }
        \item{PACE in fda:}{
          Function \code{pcaPACE} arries out a functional PCA with regularization from the 
          estimate of the covariance surface.

          Function \code{scoresPACE} estimates functional Principal Component 
          scores through Conditional Expectation (PACE).
        }
        \item{Further changes to \code{smooth.morph} and \code{landmarkreg:}}{
          \code{Smooth.morph} estimates a warping function when the target of the fit by
          registration is a functional data object.  This function has been extended 
          to work when the target for the fit and the fitted functions have different
          ranges or domains.  The warping also maps each boundary into its target
          boundary.  Simiarly \code{landmarkreg} uses a small number of discrete
          values to define the warping, and how has an extra argument, \code{x0lim},
          that defines the range of the target domain.  Since it defaults to the 
          range of the warped domain, it continues to work if not used and the
          domains have the same range.
        }
        \item{Surprisal smoothing:}{This function works with multinomial data that
          evolve over a continuum, such as the value of a latent variable in 
          psychometrics.  A multinomial observation consists of a set of 
          probabilities that are in the open interval (0,1) and sum to one.
          The surprisal value S(P_m) corresponding to a probabity P_m is 
          -log_M(P_m), where M is the number of probabities and is the base of
          the logarithm.  The inverse function is P(S_m) = M^(-S_m).
          
          Surprisal is also known as "self-information" in the field of information
          theory.  It has the characteristics of a true metric:  Surprisals can be 
          added, multiplied by positive numbers, and the difference between two 
          surprisal values mean the same thing everywhere along the information.
          continuum. The unit of the metric is called the "M-bit", the
          generalization of the familiar "bit" or "2-bit" for binary data.
          The metric property is not possessed by so-called latent  
          variables because they can be arbitrarily monotonically transformed.
          
          Smoothing surprisal data is much easier and faster than smoothing
          probabilities since surprisal values are only constrained to be 
          non-negative and are otherwise unbounded.  
          
          The function \code{smooth.surp} estimates smooth curves which fit a set of 
          surprisal values and which also satisfy the constraint that their 
          probability versions sum to one.
        }
        \item{Improvements in iterative optimisation:}{
          Many functions in the fda package optimize a fitting criterion
          iteratively.  Function \code{smooth.monotone} is an example.
          The optimisation algorithm used was a rather early design,
          and many improvements have since been made.  In most of our
          optimisations, we have switched to the algorithm to be found
          in Press, Teukolsky, Vetterling and Flannery Numerical Recipes
          volumes.  We have noticed a bit improvement in speed, are in 
          the process of upgrading all of our optimisers using this 
          approach.
        }
      }
    }
  }
}
