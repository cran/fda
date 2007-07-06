function penaltymat = fourierpen(basisobj, Lfdobj)
%FOURIERPEN computes the fourier penalty matrix for penalty LFD.
%  Arguments:
%  BASISOBJ ... a basis object
%  LFDOBJ   ... a linear differential operator object
%  Returns: 
%  PENALTYMAT ... the penalty matrix.

%  Last modified:  20 July 2006

%  check BASIS

if ~strcmp(class(basisobj), 'basis')
    error('First argument is not a basis object.')
end

%  check basis type

type = getbasistype(basisobj);
if ~strcmp(type, 'fourier')
    error ('Wrong basis type');
end

%  set up default linear differential operator

if nargin < 2, 
    Lfdobj = int2Lfd(2); 
end

%  check LFD

Lfdobj = int2Lfd(Lfdobj);

%  get highest order of derivative and check

nderiv = getnderiv(Lfdobj);
if nderiv < 0, error('NDERIV is negative');  end

%  get basis information

rangeval = getbasisrange(basisobj);
width    = rangeval(2) - rangeval(1);
params   = getbasispar(basisobj);
period   = params(1);
ratio    = round(width/period);

if width/period == ratio && isinteger(Lfdobj)

    %  Compute diagonal penalty matrix when range
    %  is a multiple of the period and LFDOBJ is integer

    penaltymat = diag(fourierpendiag(basisobj, nderiv));

else

    %  Compute penalty matrix by numerical integration in
    %  the general case.

    penaltymat = inprod(basisobj, basisobj, Lfdobj, Lfdobj);

end

%  ---------------------------------------------------------------

function pendiag = fourierpendiag(basisobj, nderiv)
%  FOURIERPENDIAG Computes fourier penalty in diagonal case

%  Last modified 2 January 2003

nbasis  = getnbasis(basisobj);
range   = getbasisrange(basisobj);
width   = range(2) - range(1);
params  = getbasispar(basisobj);
period  = params(1);
if width/period ~= round(width/period)
    error('RANGE is not a multiple of PERIOD.');
end
omega   = 2*pi/period;
halfper = period/2;
twonde  = 2*nderiv;
pendiag = zeros(1,nbasis);
if nderiv == 0
    pendiag(1) = period/2;
else
    pendiag(1) = 0;
end
j   = 2:2:nbasis-1;
fac = halfper*(j*omega/2).^twonde;
pendiag(j)   = fac;
pendiag(j+1) = fac;
pendiag = 2*pendiag./period;
pendiag = (width/period).*pendiag;


