       Updates on Functions for Functional Data Analysis in
                          R, SPLUS and Matlab

                     Jim Ramsay,  McGill University

                            6 September 2004

These changes affect only  Matlab versions of the functions :

All functions and examples have been checked for Matlab Version 7,
Release 14.  Problems in some of the Version 6.5 code were encountered,
so it is advisable to replace these functions if you have updated your
Matlab software to Version 7.

1.  The most important change is the introduction of two new objects:
   A.  The Lfd object:   This object contains information defining a
linear differential operator, and that may be nonhomogeneous.

This change in the code comes about because of work on data arising from
industrial processes where control theory can be applied.  These
processes are usually represented as nonhomogeneous linear differential
equation systems.  The new and much more powerful function for principal
differential analysis, pdacell(), was designed to estimate such systems
from data.  Once estimated, the linear differential operator can be used
to smooth data and in other ways.  The coefficients defining the linear
differential operator can each be defined separately, with their own
bases, number of basis functions, and so on.  In addition, one or more
forcing functions, each with its own weight function, may be involved.
See the preamble to function Lfd(), found in the directory @Lfd, for
more details, as well as the updated example files.

Nearly all functions are potentially affected by the introduction of the
Lfd object.  Moreover, it is now the case that few functions will work
with just an integer instead of a properly defined Lfd object.  Instead,
a new function, int2Lfd(m), converts a nonnegative integer m into a
linear differential operator object of class Lfd.

   B.  The fdPar object:  This permits us to distinguish between a
functional data object of the fd class, which has two essential slots,
the coefficient matrix and the basis, and functional parameter objects,
which are functional data objects and, in addition, also are defined by:
  -- a roughness penalty of the Lfd class,
  -- a positive real smoothing parameter, \lambda
  -- a binary indicator parameter "estimate" which determines whether
     the functional parameter is estimated (1) of held fixed (0).

The functional parameter class fdPar simplifies the argument list for
the many functions that estimate functions using a roughness penalty.
By bundling the roughness penalty and the smoothing parameter in with
the coefficient matrix and basis object, these ideas are kept together,
and the number of arguments needed in functions such as smooth_basis or
smooth_monotone is reduced.

All of the example files have been modified to take advantage of these
new objects.

These new objects make full use of Matlab's cell objects.

2.  A number of new functions are also added, including smoothing and
evaluation functions for positive, monotone, density, and warping
functions.  The registerfd() function has been reworked considerably.

3.  In the case of a B-spline basis, the function bsplinepen() can now
compute the penalty matrix exactly instead of using numerical
integration.

5.  In the case of a B-spline basis, it is now possible to have multiple
knots at a breakpoint in order to permit derivatives and even function
values to be discontinuous at chosen breakoints.  This facility has been
enabled to allow for step function inputs to dynamic systems and other
linear models.

6.  The bases that are now permitted are:
        -- constant
        -- bspline
        -- exponential
        -- fourier
        -- polygonal
        -- power

It is expected that these modifications will soon be introduced into the
S-PLUS and R code as well.

Last modified on:  6 September 2004


