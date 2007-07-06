\name{nls}
\alias{nls}
\title{
  Nonlinear Least Squares 
}
\description{
  Minor modification to 'nls' in the 'stats' package to support
  (a) estimating parameters when no 'variables' are provided, and (b)
  'data' that's a list that can not be coerced to a data.frame.  
}
\usage{
     nls(formula, data, start, control, algorithm,
         trace, subset, weights, na.action, model,
         lower, upper, ...)
}
\arguments{
  \item{formula, data, start, control, algorithm, trace, subset,
  weights, na.action, model, lower, upper, ...}{
    See help(nls, package=stats) for details.  
  }
}
\details{
  This version of 'nls' was developed on 2007.06.02 to fix some apparent
  discrepancies between documentation and behaviour.  See the examples.
  
  If the official version of 'nls' in 'stats' is modified to support the
  'fda' uses, this version might then be deleted.

  Of the options documented in the examples, the option using 'get'
  might not be supported in the future.  It may still work, however, if
  at least one argument is passed to the function via a data.frame.
  This could always be done by passing a dummy argument that is never
  used.  
}
\value{
  A list as documented in help(nls, package=stats).  
}

\seealso{
  \code{\link[stats]{nls}}
}
\examples{
# Weighted Michaelis-Menten model
# with data = a list that can not be coerced to a data.frame 
Treated <- Puromycin[Puromycin$state == "treated", ]
TreatIrreg <- with(Treated,
        list(conc1=conc[1], conc.1=conc[-1], rate=rate))
sapply(TreatIrreg, length)

# Passing arguments via 'get'
weighted.MM0 <- function(Vm, K){  
  TI <- get("TreatIrreg")
  conc <- with(TI, c(conc1, conc.1))
  resp <- TI$rate
#
  pred <- (Vm * conc)/(K + conc)
  (resp - pred) / sqrt(pred)
}

Pur.wt0 <- nls( ~ weighted.MM0(Vm, K), start = list(Vm = 200, K = 0.1),
                   trace = TRUE)
# NOTE:  'get' works but is not recommended,
# because it's too subtle and inflexible.  

# Passing arguments using a list that can not be coerced to a data.frame
weighted.MM1 <- function(resp, conc1, conc.1, Vm, K){
  conc <- c(conc1, conc.1)
#
  pred <- (Vm * conc)/(K + conc)
  (resp - pred) / sqrt(pred)
}

Pur.wt1 <- nls( ~ weighted.MM1(rate, conc1, conc.1, Vm, K), data=TreatIrreg, 
               start = list(Vm = 200, K = 0.1), trace = TRUE)

# Chambers and Hastie (1992) Statistical Models in S
# (Wadsworth, p. 537):  
# If the value of the right side [of formula] has an attribute
# called 'gradient' this should be a matrix with the number of rows
# equal to the length of the response and one column for each
# parameter.

weighted.MM.gradient <- function(resp, conc1, conc.1, Vm, K){
  conc <- c(conc1, conc.1)
#
  K.conc <- (K+conc)
  dy.dV <- conc/K.conc
  dy.dK <- (-Vm*dy.dV/K.conc)
  pred <- Vm*dy.dV
  pred.5 <- sqrt(pred)
  dev <- (resp - pred) / pred.5
  Ddev <- (-0.5*(resp+pred)/(pred.5*pred))
  attr(dev, "gradient") <- (Ddev * cbind(Vm = dy.dV, K = dy.dK))
  dev
}

Pur.wt.gradient <- nls(
   ~ weighted.MM.gradient(rate, conc1, conc.1, Vm, K), data=TreatIrreg, 
   start = list(Vm = 200, K = 0.1), trace = TRUE)

# In this example, there seems no advantage to providing the gradient.
# In other cases, there might be.  

}

% docclass is function
\keyword{models}
