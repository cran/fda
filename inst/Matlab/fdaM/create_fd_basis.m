function basisobj = create_fd_basis(fdobj)
%CREATE_FD_BASIS Creates a functional data object basis.
%  Arguments ...
%  FD ... A replicated univariate functional data object.  The replicated
%         functions are the basis functions for representing data.  A
%         popular source of such a basis would the the FD object
%         containing eigenfunctions arising from a functional PCA.
%  Returns
%  BASISOBJ  ... a functional data basis object

%  last modified 6 April 2010

%  Default basis for missing arguments: a flat line and a line with slope 1

if nargin==0
    type        = 'fd';
    rangeval    = [0,1];
    nbasis      = 2;
    params      = fd([1, 0; 0, 1],create_monomial_basis([0,1],2));
    dropind     = [];
    quadvals    = [];
    values      = {};
    basisvalues = {};

    basisobj = basis(type, rangeval, nbasis, params, ...
                     dropind, quadvals, values, basisvalues);
    return
end

if ~strcmp(class(fdobj), 'fd')
    error('Argument is not a functional data object');
end

%  extract the basis

basisobj = getbasis(fdobj);

%  extract the coefficient matrix

coefmat = getcoef(fdobj);

if length(size(coef)) > 2
    error('Argument is not a univariate functional data object');
end

%  construct basis object

type        = 'fd';
rangeval    = getbasisrange(basisobj);
nbasis      = size(coefmat,2);
params      = fdobj;
dropind     = [];
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);
