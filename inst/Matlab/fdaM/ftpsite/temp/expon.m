function basismat = expon(evalarg, ratevec, nderiv)
%EXPON Computes values of the exponentials, or their derivatives.
%  RATEVEC is a vector containing the rate constants, or mulipliers of X
%    = the exponent of e.
%  The default is the exponential function.
%  Arguments are as follows:
%  EVALARG ... A vector of values at which the polynomials are to
%              evaluated
%  RATEVEC ... a vector containing the rate constants, or mulipliers of X
%              = the exponent of e.
%  NDERIV  ... highest order derivative.  0 means only function values
%             are returned.
%  Return is a matrix with length(X) rows and NRATE columns containing
%  the values of the exponential functions or their derivatives.

%  last modified 30 January 2003

if nargin < 3, nderiv = 0;  end
if nargin < 2, ratevec = 1; end

evalarg  = evalarg(:);
n        = length(evalarg);
nrate    = length(ratevec);
basismat = zeros(n,nrate);
for irate = 1:nrate
    rate = ratevec(irate);
    basismat(:,irate) = rate.^nderiv .* exp(rate.*evalarg);
end


