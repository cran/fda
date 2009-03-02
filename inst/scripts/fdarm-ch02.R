###
###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 2.  Essential Comparisons of the Matlab and R Languages
###

library(fda)

##
## Section 2.1 A Quick Comparison of Matlab and R Syntax
##
# Names with periods:
min.fourier.basis = create.fourier.basis()

# structure of this fourier basis object:
str(min.fourier.basis)

# One component of this list of class 'basisfd'
min.fourier.basis$type

# create a vector
rng = c(0, 1)

# Access a component of a vector
rng[2]

# Logical values:
TRUE
FALSE

# T and F are variables
# by default are TRUE and FALSE, respectively,
T
F
# but can be redefined:
F = TRUE
if(F)cat('TRUE')

# alternative
F = c('Do', 'not', 'use', 'F', 'as', 'a', 'logical.')
# The following will now throw an error:
if(F)cat('TRUE')

# addition?
1 + TRUE # = 2

# Any line that is not syntactically complete
# is assumed to continue to the next line
c("like",
  'This')

"line can end in ';'";
"but not required."

# Section 2.1.2.  Using Functions
b3.4    = create.bspline.basis(norder=3, breaks=c(0, .5, 1))
fdPar3  = fdPar(b3.4, lambda=1)
fd3.4s0 = smooth.basis(0:1, 0:1, fdPar3)

class(fd3.4s0) # fdSmooth
# its 'fd' component
myfdobj. = fd3.4s0$fd

# or directly with the function call in one line
myfdobj = smooth.basis(0:1, 0:1, fdPar3)$fd

all.equal(myfdobj., myfdobj)

fd3.4s0$gcv

# specifying arguments by name not in the standard order
myfdobj = smooth.basis(y=c(1,1,2), argvals=seq(0, 1, .5), fdPar3)$fd

##
## Section 2.2 Singleton Index Issues
##
temp = matrix(c(1,2,3,4),2,2)
class(temp)
class(temp[,1])
temp[,1] # not a matrix
temp[,1, drop=FALSE] # still a matrix

index = 1:2
temp[, index] # matrix
index = 1
temp[, index] # not a matrix

index = 1:2
temp[, index][, 1] # OK
index = 1
temp[, index][, 1] # Error

# 3-d array
A      = array(1:2, dim=c(1, 2, 1))
a1     = A[, 1, ] # scalar
dim(a1)= dim(A)[-2]
a1 # 1 x 1 matrix

##
## Section 2.3 Classes and Objects in R and Matlab
##
default.fd = fd()
is.list(default.fd)

class(default.fd$basis)

plot(default.fd)
# same as
plot.fd(default.fd)

attributes(default.fd)

coef(default.fd)
# same as
default.fd$coefs

# fdPar example
rangeval  = c(-3,3)
x         =  rnorm(50)
x[x < -3] = -2.99
x[x >  3] =  2.99
basisobj  = create.bspline.basis(rangeval, 11)
Wfd0      = fd(matrix(0,11,1), basisobj)
WfdParobj = fdPar(Wfd0)

coef(WfdParobj)
WfdParobj$coefs # NULL
WfdParobj$fd$coefs # OK

# str = structure
str(WfdParobj)

# All methods for the generic function 'coef' in attached packages
methods(coef)
library('nlme')
# add 'coef' methods for objects defined in the 'nlme' package
methods(coef)

methods(class='fd') # methods for objects of class 'fd'

##
## Section 2.4 More to Read
##
# should open a locally install page in a browser
help.start()


