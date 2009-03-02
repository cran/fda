function  coef = project_basis(y, argvals, basisobj, penalize)
% PROJECT_BASIS Project discrete values on to a set of basis fns.
%  Arguments for this function:
%
%  Y        ... an array containing values of curves
%               If the array is a matrix, rows must correspond to argument
%               values and columns to replications, and it will be assumed
%               that there is only one variable per observation.
%               If Y is a three-dimensional array, the first dimension
%               corresponds to argument values, the second to replications,
%               and the third to variables within replications.
%               If Y is a vector, only one replicate and variable are assumed.
%  ARGVALS  ... A vector of argument values.  This must be of length
%               length(Y) if Y is a vector or size(Y)(1) otherwise.
%  BASISOBJ  ... A basis.fd object
%  PENALIZE  ... If nonzero, a penalty matrix is used to deal with
%                  a singular basis matrix.  Normally this is not needed.
%
%  Returns a coefficient vector or array. The first dimension is the number
%     of basis functions and the other dimensions (if any) match
%     the other dimensions of Y.
%

%  Last modified:  20 July 2006

if nargin < 4
    penalize = 0;
end

%  Calculate the basis and penalty matrices, using the default
%   for the number of derivatives = the penalty.
basismat = getbasismatrix(argvals, basisobj);
Bmat     = basismat' * basismat;
nbasis   = getnbasis(basisobj);
if penalize ~= 0
    %  Add a very small multiple of the identity to penmat
    %   and find a regularization parameter
    penmat = eval_penalty(basisobj, 0);
    penmat = penmat + 1e-10 .* max(max(penmat)) .* eye(nbasis);
    lambda1 = (0.0001 .* sum(diag(Bmat)))./sum(diag(penmat));
    Cmat = Bmat + lambda1 .* penmat;
else
    Cmat = Bmat;
end
%  Do the fitting by a simple solution of the
%    equations taking into account smoothing
ydim = size(y);
if(length(ydim) <= 2)
    Dmat = basismat' * y;
    coef = symsolve(Cmat,Dmat);
else
    nvar = ydim(3);
    coef = zeros([nbasis, ydim(2), nvar]);
    for ivar = 1:nvar
        Dmat = basismat' * y(:,:,ivar);
        coef(:,:,ivar) = symsolve(Cmat,Dmat);
    end
end

