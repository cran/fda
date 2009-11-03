function bibasisobj =  create_product_basis(sbasisobj, tbasisobj)
%CREATE_PRODUCT_BASIS creates a basis object for bi-variate functions.
%  A product basis basis object is constructed from
%  two univariate bases for expanding a bivariate function.  
%  Each univariate basis is contained in a 'basis' object.  
%  The basis for the bivariate functional data object consists
%  of all possible pairs of the two sets of univariate basis functions.
%  Arguments
%  SBASISOBJ ... a functional data basis object for the first  argument s
%  TBASISOBJ ... a functional data basis object for the second argument t

%  Returns
%  BASISOBJ ... a functional data object

%  last modified 20 July 2006

if nargin == 0
    bibasisobj.sbasisobj = basis('bspline',[0,1],3,0.5);
    bibasisobj.tbasisobj = basis('bspline',[0,1],3,0.5);
    bibasisobj = class(bibasisobj, 'basis');
    return;
end

if ~isa_basis(sbasisobj)
    error('Argument SBASISOBJ must be of basis class');
end
if ~isa_basis(tbasisobj)
    error('Argument TBASISOBJ must be of basis class');
end

type      = 'product';
rangevals = getbasisrange(sbasisobj);
rangevalt = getbasisrange(tbasisobj);
nbasiss   = getnbasis(sbasisobj);
nbasist   = getnbasis(tbasisobj);
nbasis    = nbasiss*nbasist;
params.sbasis = sbasisobj;
params.tbasis = tbasisobj;

bibasisobj = bibasis(type, rangevals, rangevalt, nbasis, params);


