se###
###
### Ramsey & Silverman (2006) Functional Data Analysis, 2nd ed. (Springer)
###
### ch. 19.  Fitting Differential Equations to Functional Data:
###          Principal Differential Analysis 
###
library(fda)
# p. 327 
##
## Section 19.1.  Introduction
##


##
## Section 17.2.  Oil refinery data 
##
#. p. 298, Figure 17.1.  Reflux flow and Tray 47 level in a refinery 

str(refinery)
sapply(refinery, class)

with(growth, matplot(age, hgtf[, sel], type="b"))

refOrder <- 4
(refKnots <- c(0, 67/2, 67, 67, 67, 67+(1:4)*(193-67)/4))
(ref.nKnots <- length(refKnots))
# p. 299:  "These knot choices imply a total of eleven basis functions."  
(ref.nBasis <- ref.nKnots+refOrder-2)

# I should work fda-ch14 before continuing here.  
                   
