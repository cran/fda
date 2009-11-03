function basisobj = create_monom_basis(rangeval, nbasis, exponents)
%  CREATE_MONOM_BASIS  Creates a monomial basis:, x^i_1, x^i_2, ...
%  The exponents in this version must be nonnegative integers
%  Argument:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a 
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS    ... number of basis functions
%  EXPONENTS ... an array of NBASIS nonnegative integer exponents
%                by default this is 0:(NBASIS-1)
%  Return:
%  BASIS  ... a functional data basis object of type 'monom'

%  last modified 24 October 2003

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not positive.');
    end
    rangeval = [0,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

if nargin < 2, nbasis = 2;    end
if nargin < 3, exponents = 0:(nbasis-1); end

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
type   = 'monom';
params = exponents;

basisobj = basis(type, rangeval, nbasis, params);

