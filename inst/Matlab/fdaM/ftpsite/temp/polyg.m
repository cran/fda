function basismat = polyg(evalarg, argvals, nderiv)
%  POLYG Evaluates the basis for a linear interpolant or its first derivative.
%  It calls function spcol.
%  Arguments are as follows:
%  EVALARG ... A vector of values at which the spline functions are to
%              evaluated
%  ARGVAL  ... a STRICTLY INCREASING sequence of argument values.
%  NDERIV  ... Either 0 or 1.  0 means only function values
%  Return is a matrix with length(EVALARG) rows and number of columns equal to
%             number of argument values

%  last modified 20 July 2006

if nargin < 3
    nderiv = 0;
end

argvals  = argvals(:);
nargvals = length(argvals);
range    = [argvals(1), argvals(nargvals)];

evalarg = evalarg(:);
if (max(evalarg) > max(argvals)) || (min(evalarg) < min(argvals)) 
    error('ARGVALS do not span the values of EVALARG.');
end

if (min(diff(argvals)) <= 0 )
    error('Break-points are not strictly increasing');
end

if (~(nderiv == 0 || nderiv == 1))
    error('NDERIV is neither 0 nor 1.');
end

nbasis = length(argvals);

basis = create_bspline_basis(range, nbasis, 2, argvals);

basismat = eval_basis(evalarg, basis, nderiv);

