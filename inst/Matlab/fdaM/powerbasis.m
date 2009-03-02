function powermat = powerbasis(evalarg, exponents, nderiv)
%POWERBASIS computes values of monomials, or their derivatives.
%  The powers of EVALARG are the NBASIS nonnegative integers in EXPONENTS.
%  The default is 1, meaning EVALARG itself.
%  Arguments are as follows:
%  EVALARG   ... A vector of values at which the polynomials are to
%                evaluated.
%  EXPONENTS ... vector of exponents
%  NDERIV    ... order of derivative to be returned.
%  Return is:
%  A matrix with length(EVALARG) rows and NBASIS columns containing
%    the values of the monomials or their derivatives

%  last modified 20 July 2006

evalargdim = size(evalarg);
if evalargdim(1) > 1 && evalargdim(2) > 1
    error('Argument EVALARG is not a vector.');
end
evalarg = evalarg(:);
n = length(evalarg);

% set default arguments

if nargin < 3, nderiv = 0; end

nbasis = length(exponents);

powermat = zeros(n,nbasis);
if nderiv == 0
    for ibasis=1:nbasis
        powermat(:,ibasis) = evalarg.^exponents(ibasis);
    end
else
    if any(exponents - nderiv < 0) && any(evalarg == 0)
        error('A negative exponent is needed and an argument value is 0.');
    else
        for ibasis=1:nbasis
            degree = exponents(ibasis);
            fac = degree;
            for ideriv=2:nderiv
                if degree == ideriv - 1
                    fac = 0;
                else
                    fac = fac*(degree-ideriv+1);
                end
            end
            powermat(:,ibasis) = fac.*evalarg.^(degree-nderiv);
        end
    end
end

