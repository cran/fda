function bibasisobj = bibasis(basistype, srangeval, trangeval, ...
    nbasis, params)
%BIBASIS  Creates a tensor product functional data basis over a rectangle.
%  Arguments
%  BASISTYPE ...  a string indicating the type of basis.  
%                 This may be one of:
%                'product', 'prod', 'Product', 'Prod',
%                'history', 'hist', 'History', 'Hist', 'Historical',
%  SRANGEVAL ... an array of length 2 containing the lower and upper
%                boundaries for the rangeval of values for the
%                first argument "s"
%  TRANGEVAL ... an array of length 2 containing the lower and upper
%                boundaries for the rangeval of values for the
%                first argument "t"
%  NBASIS   ... the number of basis functions
%  PARAMS   ... If the basis is 'product', this is a struct object with
%               two slots:  sbasis and tbasis.  
%               PARAMS.SBASIS contains a basis object for argument "s"
%               PARAMS.TBASIS contains a basis object for argument "t"
%               If the basis is 'history', PARAMS is a vector created
%               by the command
%               params = [reshape(eleNodes,nbasis*3,1); Si; Ti].
%               The first NBASIS*3 elements are indices of nodes that
%               are the vertices of the triangular elements, and
%               the remainder are "s" and "t" values for the nodes.
%  Returns
%  BIBASIS_fd  ... a tensor product basis object

%  Specific types of bases may be set up more conveniently using functions
%  CREATE_PRODUCT_BASIS  ...  creates a tensor product basis

%  last modified 7 March 2011

if nargin==0
    bibasisobj.type      = 'product';
    bibasisobj.srangeval = [0,1];
    bibasisobj.trangeval = [0,1];
    bibasisobj.nbasis    = 1;
    bibasisobj.params.sbasis = create_constant_basis([0,1]);
    bibasisobj.params.tbasis = create_constant_basis([0,1]);
    bibasisobj = class(bibasisobj, 'bibasis');
    return;
end

if isa(basistype, 'bibasis')
    bibasisobj = basistype;
    return;
end

%  check PARAMS vector according to the basis type
switch basistype
    case 'product'
        if ~strcmp(params, 'struct')
            error ('PARAMS is not a struct object.');
        end
        if ~strcmp(params.sbasis, 'basis')
            error('PARAMS.SBASIS is not a basis object.');
        end
        if ~strcmp(params.tbasis, 'basis')
            error('PARAMS.TBASIS is not a basis object.');
        end
    case 'history'
        if ~strcmp(params, 'double')
            error('PARAMS is not a vector.');
        end
        nparams = length(params);
        if nparams <= nbasis
            error('PARAMS is not longer than NBASIS.');
        end
        ncoords = nparams - nbasis;
        if 2*floor(ncoords/2) ~= 2
            error('Length of PARAMS after NBASIS not even.');
        end
    otherwise
        error('Unrecognizable bivariate basis');
end

bibasisobj.type      = basistype;
bibasisobj.srangeval = srangeval;
bibasisobj.trangeval = trangeval;
bibasisobj.nbasis    = nbasis;
bibasisobj.params    = params;

bibasisobj = class(bibasisobj, 'bibasis');

