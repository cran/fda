###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 3.  How to specify basis systems for building functions
###
library(fda)
##
## Section 3.1 Basis Function Systems for Constructing Functions
##
const.basis   = create.constant.basis(rangeval=0:1)
monom.basis   = create.monomial.basis(rangeval=0:1, nbasis=1)
fourier.basis = create.fourier.basis(rangeval=0:1,
                      nbasis=5, period=1)
bspline.basis = create.bspline.basis(rangeval=0:1,
             nbasis=5, norder=2, breaks=seq(0, 1, length=5) )

##
## Section 3.2 Fourier Series for Periodic Data and Functions
##
daybasis65 = create.fourier.basis(c(0,365), 65)

T.         = 500
daybasis.T = create.fourier.basis(c(0,365), 3, period=T.)

zerobasis  = create.fourier.basis(c(0, 365), 65, dropind=1)
zerobasis. = daybasis65[2:65]
all.equal(zerobasis, zerobasis.)

fourier.basis3 = create.fourier.basis(rangeval=c(0, 1), nbasis=3)

help(create.fourier.basis)

##
## Section 3.3 Spline series for Non-periodic Data and Functions
##

#  section 3.3.3 Examples
bspline4 = create.bspline.basis(breaks=c(0, .5, 1))
knots(bspline4, interior=FALSE)

bspline2 = create.bspline.basis(breaks=c(0, .5, 1), norder=2)
knots(bspline2, interior=FALSE)

bspline2.2 = create.bspline.basis(breaks=c(0, .5, .5, 1), norder=2)
knots(bspline2.2, interior=FALSE)

bspline4.3 = create.bspline.basis(breaks=c(0, .5, .5, .5, 1))
knots(bspline4.3, interior=FALSE)

#  section 3.3.4 B-Splines
splinebasis = create.bspline.basis(c(0,10), 13)

# Figure 3.1
plot(splinebasis, xlab='t', ylab='Bspline basis functions B(t)', las=1)

# Figure 3.2
basis2 = create.bspline.basis(c(0,2*pi), 5, 2)
basis3 = create.bspline.basis(c(0,2*pi), 6, 3)
basis4 = create.bspline.basis(c(0,2*pi), 7, 4)

theta    = seq(0, 2*pi, length=201)
sin.theta= sin(theta)

sin2 = Data2fd(theta, sin.theta, basis2)
sin3 = Data2fd(theta, sin.theta, basis3)
sin4 = Data2fd(theta, sin.theta, basis4)

sin2.theta = predict(sin2, theta)
sin3.theta = predict(sin3, theta)
sin4.theta = predict(sin4, theta)

sinRng = range(sin2.theta)
pi3    = ((1:3)*pi/2)

op = par(mfrow=c(3,2), mar=c(3,4,2,2)+.1)

plot(theta, sin2.theta, type='l', ylim=sinRng, xlab='', ylab='Order = 2',
     main='sine(t)' )
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin2.theta = predict(sin2, theta, 1)
plot(theta, Dsin2.theta, type='l', ylim=sinRng, xlab='', ylab='',
     main='D sine(t)')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

plot(theta, sin3.theta, type='l', ylim=sinRng, xlab='', ylab='Order = 3')
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin3.theta = predict(sin3, theta, 1)
plot(theta, Dsin3.theta, type='l', ylim=sinRng, xlab='', ylab='')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

plot(theta, sin4.theta, type='l', ylim=sinRng, xlab='t', ylab='Order = 4')
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin4.theta <- predict(sin4, theta, 1)
plot(theta, Dsin4.theta, type='l', ylim=sinRng, xlab='t', ylab='')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

par(op)

splinebasis6 = create.bspline.basis(c(0,10), 15, 6)

# Knots not equally spaced
spline.unequal.knots = create.bspline.basis(breaks=c(0, .7, 1))
knots(spline.unequal.knots, interior=FALSE)

plot(splinebasis6[7])
plot(splinebasis6[1:2])

t10 = seq(0, 10, .1)
spline6.10 = predict(splinebasis6, t10)
plot(t10, rowSums(spline6.10))

heightbasis  = with(growth, create.bspline.basis(c(1, 18), 35,6,age))
# or
heightbasis. = create.bspline.basis(norder=6, breaks=growth$age)
all.equal(heightbasis, heightbasis.)

help(create.bspline.basis)

##
## Section 3.4 Constant, Monomial and Other Bases
##
conbasis = create.constant.basis(c(0,1))
conb     = create.constant.basis()
all.equal(conbasis, conb)

monbasis = create.monomial.basis(c(0,1), 4)
monb     = create.monomial.basis(nbasis=4)
all.equal(monbasis, monb)

##
## Section 3.5 Methods for Functional Basis Objects
##
methods(class='basisfd')

print(monb)
print.basisfd(monb) # the same

summary(monb)

monbasis==monb

is(monb, 'basisfd')
inherits(monb, 'basisfd')

monb$params
monb$params = 2*(0:3)

t.1       = seq(0, 1, .1)
monbmat.1 = eval.basis(t.1, monb)
Dmonbmat.1= eval.basis(t.1, monb, 1)

monbmat.1. = predict(monb, t.1)
all.equal(monbmat.1, monbmat.1.)

Dmonbmat.1.= predict(monb, t.1, 1)
all.equal(Dmonbmat.1, Dmonbmat.1.)

##
## Section 3.6 The Structure of the basisfd or basis Class
##

help(basisfd)

##
## Section 3.7 Some Things to Try
##
# (exercises for the reader)
