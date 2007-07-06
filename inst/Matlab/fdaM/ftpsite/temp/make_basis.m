function basisobj = make_basis(rangeval, nresol, periodic, nderiv)
%MAKE_BASIS sets up a simple basis.
%  If argument PERIODIC is nonzero, the basis is of type 'fourier',
%  else it is of type 'bspline'.
%  Argument RANGEVAL determines the range of argument values.  If
%  the basis type is 'fourier', this also determines the period.
%  The number of basis functions is jointly determined by arguments
%  NRESOL and NDERIV.
%  For a 'fourier' basis, NBASIS = NRESOL.
%  For a 'bspline' basis, NBASIS = NRESOL + NDERIV + 4
%
%  Arguments are as follows:
%
%  RANGEVAL ... A vector of length 2 giving the lower and upper limits
%               on argument values, respectively.
%  NRESOL   ... The resolution required in the functions.  This means
%               the maximum number of features or events that could in
%               principle be represented in each observation.  Features
%               are things like peaks, valleys, zero crossing, or plateaus.
%               The width of a single feature may naturally not be less than 
%               the minimum difference between two successive argument values.
%               Roughly speaking, NRESOL is the width of the interval in 
%               RANGEVAL divided by the width of the narrowest feature that
%               requires representation.  NRESOL cannot logically exceed
%               the number N of argument values, and for noisy data it should
%               be considerably less. A reasonable rule of thumb is to
%               NRESOL to 1/3 of the number of sampling points, if these
%               are more or less evenly spaced.
%  PERIODIC ... Is 0 if the functions are not periodic, and nonzero otherwise.
%               It is set to 0 by default.
%  NDERIV   ... The highest order of derivative that is needed for the 
%               functional data object.  This is set to 2 by default.
%
%  Another option for a simple basis may be the polygonal basis of type 'polyg'
%    made by function CREATE_POLYGONAL_BASIS.  Only use this option if no 
%    derivatives will be required.
%
%  Last modified  1 December 2000

if nargin < 4
   nderiv = 2;
end
if nargin < 3
   periodic = 0;
end
if nargin < 2
   error('Argument NRESOL must be supplied.');
end
nresol = floor(nresol);
if nresol < 0
   error('Argument NRESOL must be nonnegative.');
end
if length(rangeval) ~= 2
   error('Argument RANGEVAL must be of length 2.');
end
width = rangeval(2) - rangeval(1);
if  width <= 0
   error('Values in argument RANGEVAL must be strictly increasing.');
end

if periodic
   %  periodic basis
   nbasis   = nresol;
   basisobj = create_fourier_basis(rangeval, nbasis, width);
else
   %  B-spline basis
   nbasis   = nresol + nderiv + 4;
   basisobj = create_bspline_basis(rangeval, nbasis);
end
