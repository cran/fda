function monommat = monom(evalarg, exponents, nderiv)
%  MONOM  Values of monomials, or their derivatives.
%  The powers of EVALARG are the NBASIS nonnegative integers in EXPONENTS.
%  The default is 1, meaning EVALARG itself.
%  Arguments are as follows:
%  EVALARG   ... array of values at which the polynomials are to
%                evaluated
%  EXPONENTS ... array of nonnegative integer exponents of EVALARG
%  NDERIV    ... order of derivative to be returned.
%  Return is:
%  A matrix with length(EVALARG) rows and NBASIS columns containing
%    the values of the monomials or their derivatives

%  last modified 30 January 2003

% set default arguments

if nargin < 3
    nderiv = 0;
end

if nargin < 2
    exponents = 1;
end

evalarg = evalarg(:);
n = length(evalarg);

nbasis = length(exponents);

%  check whether exponents are nonnegative integers

for ibasis=1:nbasis
    if exponents(ibasis) - round(exponents(ibasis)) ~= 0
        error('An exponent is not an integer.');
    end
    if exponents(ibasis) < 0
        error('An exponent is negative.');
    end
end

% check if there are duplicate exponents

if min(diff(sort(exponents))) <= 0
    error('There are duplicate exponents.');
end

monommat = zeros(n,nbasis);

if nderiv == 0
    %  use the recursion formula to compute monomnomial values
    for ibasis=1:nbasis, monommat(:,ibasis) = evalarg.^exponents(ibasis); end
else
    for ibasis=1:nbasis
        degree = exponents(ibasis);
        if nderiv <= degree 
            fac = degree;
            for ideriv=2:nderiv
                fac = fac*(degree-ideriv+1);
            end
            monommat(:,ibasis) = fac.*evalarg.^(degree-nderiv);
        end
    end
end

