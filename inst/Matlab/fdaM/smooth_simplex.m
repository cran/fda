function [Wfd, Warray, Parray] = ...
   smooth_simplex(argvals, y, WfdPar, conv, iterlim, dbglev)
%SMOOTH_SIMPLEX smooths the relationship of observations 
%  in a d-dimensional simplex distributed ARGVALS.   The data are in
%  an N by NCURVES by D array, where N is the length of ARGVALS,
%  NCURVES is the number of replications, and D is the dimension of the
%  simplex.  For example, if observations are in an equilateral triangle,
%  then D = 3.
%  An error sum of squares criterion is minimized by the Gauss-Newton
%  method used in Matlab function lsqnonlin.
%
%  The fitting criterion is penalized least squares.
%
%  Arguments are:
%  ARGVALS ... A vector of argument values
%  Y       ... An N by NCURVES by D array of data
%  WFDPAR  ... functional parameter object defining initial smoothing
%              functions.  This is a D-variate functional data object.
%  ITERLIM ... iteration limit for scoring iterations
%  CONV    ... convergence criterion
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFD      ... functional data basis object for the W-functions w_v(t).
%  WARRAY   ... N by NCURVES by D array of values of the W-functions
%               at the argument values in ARGVALS
%  PARRAY   ... N by NCURVES by D array of values of the P-functions
%               at the argument values in ARGVALS

%  Last modified 30 June 2010

if nargin < 2
   error('Less than two arguments are supplied.');
end

%  check ARGVALS

[argvals, n] = argcheck(argvals);

%  check Y

[y, ncurve, d, ndim] = ycheck(y, n);

if ndim < 3
    error('Data array Y does not have three dimensions.');
end

%  get basis information

if isa_fd(WfdPar) || isa_basis(WfdPar)
    WfdPar = fdPar(Wfd, 2, 0);
end

Wfd   = getfd(WfdPar);
Wcoef = getcoef(Wfd);
[nbasis,N,dim] = size(Wcoef);
lambda = getlambda(WfdPar);

if N ~= ncurve
    error('2nd dim. of Wfd coef. matrix not equal to NCURVE.');
end
if dim ~= d
    error('3rd dim. of Wfd coef. matrix not equal to D.');
end

basisobj = getbasis(Wfd);

phimat = eval_basis(argvals, basisobj);
zmat   = zerobasis(d);

%  set some default arguments and constants

if nargin < 6, dbglev  = 1;    end
if nargin < 5, iterlim = 50;   end
if nargin < 4, conv    = 1e-4; end

%  check LFDOBJ

Lfdobj = getLfd(WfdPar);
Lfdobj = int2Lfd(Lfdobj);

dbgwrd  = dbglev > 1;

%  initialize matrix Kmat defining penalty term

if lambda > 0
  Kmat = eval_penalty(basisobj, Lfdobj);
  [V,D] = eig(full(Kmat));
  D = diag(D);
  [Dsort, Isort] = sort(D);
  Vsort = V(:,Isort);
  nderiv = getnderiv(Lfdobj);
  ind = (nderiv+1):nbasis;
  Vsort = Vsort(:,ind);
  Dsort = Dsort(ind);
  Kfac = Vsort*diag(sqrt(Dsort));
  Kfac = sqrt(lambda).*Kfac;
end

%  options for nonlinear LS minimization

if dbglev > 0
    lsOptions = optimset('Display','iter', 'Jacobian','on', ...
                         'MaxIter', iterlim);
else
    lsOptions = optimset('Display','off', 'Jacobian','on', ...
                         'MaxIter', iterlim);
end    

%  --------------------------------------------------------------------
%              loop through and curves
%  --------------------------------------------------------------------

carray = zeros(nbasis,ncurve,d);
Parray = zeros(n,ncurve,d);
Warray = zeros(n,ncurve,d);
for icurve=1:ncurve
    ymat = squeeze(y(:,icurve,:));
    yfd  = Wfd(icurve,:);
    cmat = squeeze(getcoef(yfd));
    bmat = cmat*zmat;
    bvec0 = reshape(bmat,nbasis*(d-1),1);
    [resvec, DSSE] = SSEfnpos(bvec0, n, nbasis, d, ymat, phimat, zmat, ...
                                   Kfac, lambda);
    bvec  = lsqnonlin(@SSEfnpos, bvec0, [], [], lsOptions, ...
                      n, nbasis, d, ymat, phimat, zmat, Kfac, lambda);
    [resvec, DSSE] = SSEfnpos(bvec, n, nbasis, d, ymat, phimat, zmat, ...
                                   Kfac, lambda);
    bmat = reshape(bvec,nbasis,d-1);
    cmat = bmat*zmat';
    carray(:,icurve,:) = cmat;
    Wfdi = fd(cmat,basisobj);
    Wmat = eval_fd(argvals, Wfdi);
    Emat = exp(Wmat);
    Svec = sum(Emat,2);
    Pmat = Emat./(Svec*ones(1,d));
    Parray(:,icurve,:) = Pmat;
    Warray(:,icurve,:) = Wmat;
end

Wfd = fd(carray, basisobj);





