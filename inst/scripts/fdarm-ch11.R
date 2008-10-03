###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 11  Functional Models and Dynamics
###
library(fda)

##
## Section 11.1  Introduction to Dynamics
##
#  Figure 11.1





##
## Section 11.2 Principal Differential Analysis for Linear Dynamics
##
#  (no computations in this section)

##
## Section 11.3 Principal Differential Analysis of the Lip Data
##
# Figure 11.2
bwtlist = list(fdPar(lipbasis,2,0),
    fdPar(lipbasis,2,0))
pdaList = pda.fd(lipfd,bwtlist)

plot.pda(pdaList)
dfd <- 0.25*pdaList$bwtlist[[2]]$fd^2
             - pdaList$bwtlist[[1]]$fd
dfd$fdnames = list(’time’,’rep’,’discriminant’)

# Figure 11.3 ????



# Figure 11.4
pda.overlay(pdaList)



##
## Section 11.4 PDA of the Handwriting Data
##
xfdlist = list(fdafd[,1],fdafd[,2],fdafd[,3])

pdaPar = fdPar(fdabasis,2,0)
pdaParlist = list(pdaPar, pdaPar)
bwtlist = list( list(pdaParlist,pdaParlist),
    list(pdaParlist,pdaParlist) )
pdaList = pda.fd(xfdlist, bwtlist)


eigen.pda(pdaList)

# Figure 11.5









##
## Section 11.5 Registration and PDA
##
lipreglist = landmarkreg(lipfd, as.matrix(lipmarks),
    lipmeanmarks, WfdPar)
Dlipregfd = register.newfd(deriv.fd(lipfd,1),
    lipreglist$warpfd)
D2lipregfd = register.newfd(deriv.fd(lipfd,2),
    lipreglist$warpfd)
xfdlist = list(-Dlipregfd,-lipreglist$regfd)
lipregpda = fRegress( D2lipregfd, xfdlist, bwtlist)



##
## Section 11.6 Details for pda.fd, eigen.fd, pda.overlay
##              and register.newfd
##



##
## Section 11.7 Some Things to Try
##
# (exercises for the reader)

##
## Section 11.8  More to Read
##
