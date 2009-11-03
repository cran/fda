function powerfd = mpower(fdobj, m)
%  A positive integer power of a functional data object with a B-spline
%  basis.
%  Power M is tested for being a positive integer, and, if not, 
%  an error message is issued.  
%  The basis is tested for being a B-spline basis.
%  If M and the basis of FDOBj pass these tests, then  the function
%  sets up a new spline basis with the same knots but with an order
%  that is M-1 higher than the basis for FDOBJ, so that the 
%  order of differentiability for the new basis is appropriate.
%  The power of the values of the function over a fine mesh are computed,
%  and these are fit using the new basis.
%
%  Powers should be requested with caution, however.  If there is
%  strong local curvature in FDOBJ, and if its basis is just barely
%  adequate to capture this curvature, then this function may introduce
%  considerable error into the values of the power over this local area.
%
%  If a power of a functional data object is required for which the
%  basis is not a spline, it is better to either re-represent the
%  function in a spline basis, or, perhaps even better, to do the 
%  math required to get the right basis and interpolate function
%  values over a suitable mesh.  This is especially true for fourier
%  bases.

%  Last modified 26 February 2009

if nargin < 2,  error('Number of arguments is less than two.');  end

%  Test M for being a positive integer

if ~isinteger(m) || m < 0
    error('M is not a positive integer.');
end

basisobj = getbasis(fdobj);

%  test the basis for being of B-spline type

if ~strcmp('bspline',getbasistype(basisobj))
    error('FDOBJ does not have a spline basis.');
end

nbasis = getnbasis(basisobj);
rangeval = getbasisrange(basisobj);
interiorknots = getparams(basisobj);
norder = nbasis - length(interiorknots);

coefmat = getcoef(fdobj);
coefd = size(coefmat);
    ncurve = coefd(2);
if length(coefd) == 2
    nvar = 1;
else
    nvar = coefd(3);
end
    
%  M == 0:  set up a constant basis and return the unit function(s)

if m == 0
    newbasis = create_constant_basis(rangeval);
    if nvar == 1
      powerfd = fd(ones(1,ncurve), newbasis);
    else
      powerfd = fd(ones(1,ncurve,nvar), newbasis);
    end
    return;
end

%  M == 1:  return the function

if m == 1
    powerfd = fdobj;
    return;
end

%  M > 1:

breaks = [rangeval(1), interiorknots, rangeval(2)];
newbasis = create_bspline_basis(rangeval, nbasis+m-1, norder+m-1, breaks);
tval = linspace(rangeval(1),rangeval(2),10*nbasis+1);
fmat = eval_fd(tval, fdobj);
newfdPar = fdPar(newbasis, max([nbasis-2,0]),1e-10);
powerfd  = smooth_basis(tval, fmat, newfdPar);

