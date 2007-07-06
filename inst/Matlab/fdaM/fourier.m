function basismat = fourier(evalarg, nbasis, period, nderiv)
%  FOURIER  Computes the NDERIV derivative of the Fourier series basis
%    for NBASIS functions with period PERIOD, these being evaluated
%    at values in vector EVALARG.
%  Returns an N by NBASIS matrix BASISMAT of function values

%  last modified 30 July 2003

evalarg = evalarg(:); %  ensure that EVALARG is a column vector
n       = length(evalarg);
onen    = ones(n,1);
range   = [min(evalarg),max(evalarg)];

%  set default number of basis functions

if nargin < 2
    nbasis = n;
end

%  set default period

if nargin < 3
    period = range(2) - range(1);
end

%  set default order of derivative

if nargin < 4
    nderiv = 0;
end

%  check argument values

if nbasis <= 0,  error('NBASIS not positive');  end
if period <= 0,  error('PERIOD not positive');  end
if nderiv <  0,  error('NDERIV negative');      end

%  make number of basis functions odd if required

if 2*floor(nbasis/2) == nbasis
    nbasis = nbasis + 1;
end

%  set up the basis matrix

basismat = zeros(n,nbasis);

%  set up some constants

omega = 2*pi/period;
fac   = omega .* evalarg;

if nderiv == 0
    %  The fourier series itself is required.
    basismat(:,1) = 0.7071068;
    j    = 2:2:(nbasis-1);
    args = fac * (j./2);
    basismat(:,j)   = sin(args);
    basismat(:,j+1) = cos(args);
else
    %  A derivative of the fourier series is required.
    basismat(:,1) = 0.0;
    if nderiv == floor(nderiv/2)*2
        mval  = nderiv/2;
        ncase = 1;
    else
        mval  = (nderiv-1)/2;
        ncase = 2;
    end
    j    = 2:2:(nbasis-1);
    fac  = onen * (((-1).^mval).*((j./2).*omega).^nderiv);
    args = (omega.*evalarg) * (j./2);
    if ncase == 1
        basismat(:,j)   =  fac .* sin(args);
        basismat(:,j+1) =  fac .* cos(args);
    else
        basismat(:,j)   =  fac .* cos(args);
        basismat(:,j+1) = -fac .* sin(args);
    end
end
basismat = basismat./sqrt(period./2);
