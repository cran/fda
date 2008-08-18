as.function.basisfd <- function (x, variable='x', ...){
##
## 1.  Create an environment with the basisfd 
##
  fenv <- new.env()
#  
  assign('basisobj', x, envir=fenv)  
##
## 1.  Create an empty function skeleton
##  
  f.ch <- paste('function(', variable,
                ', ...){eval.basis(', variable,
                ', basisobj, Lfdobj=0)}', sep='')
  fp <- parse(text=f.ch)

  
  fc <- call('f', f.ch)
#
# object of class function but without a body  
  f0 <- eval(fp, fenv) 
# f0 = empty function   
##
## 2.  Create the body
##
  

  
#  fCh <- paste('eval.basis(', variable, ', ', x, ')')
#  body(f0) <- parse(text=fCh)
#  f0
  'comiing soon' 
}

as.function.fd <- function (x, variable='x', ...){
  'coming soon' 
}

as.function.fdSmooth <- function (x, variable='x', ...){
  'coming soon' 
}

