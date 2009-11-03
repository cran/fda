function [smooth_fd, laplace_fd] = ....
                  smooth_FEM_fd(data, fdobj, lambda, desmat)
% SMOOTH_FEM_FD Smooth data observed at a subset of nodes, possibly
%  using covariate values in DESMAT as additional contributors to
%  the fit to the data values.
%
%     Arguments:
% FOBJ   ... a functional data object having a basis of the FEM type.
% DATA   ... a n-by-2 new set of observations.  DATA(:,1)
%            indexes the points (FELSPLOBJ.POINTS) at which the 
%            values in DATA(:,2) are observed.  If DATA is n by 1,
%            the data are assumed to be associated with the first n
%            nodes.
% LAMBDA ... a scalar smoothing parameter, MUST be positive
% DESMAT ... a design matrix
%
%     Output:
% SMOOTH_FD  ...  A FD object of the FEM type defined by the coefficient
%                 vector resulting from smoothing
% LAPLACE_FD ...  A FD object of the FEM type for the value of the 
%                 Laplace operator
%

% Last modified on 16 June 2011 by Jim Ramsay

%  assign defaults to missing arguments

if nargin < 4,  desmat = [];     end
if nargin < 3,  lambda = 1e-12;  end

if nargin < 2
    error('Only one argument supplied.  At two required.');
end

%  ---------------------------------------------------------------
%  check arguments
%  ---------------------------------------------------------------

%  FDOBJ

if ~isa_fd(fdobj) && ~isa_basis(fdobj)
   error('FDOBJ is not a FD object and not an FEM basis object.');
end

%  BASISOBJ defining FDOBJ

if isa_fd(fdobj)
    basisobj = getbasis(fdobj);
else
    basisobj = fdobj;
end

basistype = getbasistype(basisobj);
if ~strcmp(basistype, 'FEM')
    error('FDOBJ does not have a basis of the FEM type.');
end

%  LAMBDA

if  ~isa(lambda,'double')
   error('LAMBDA is not numeric')
elseif size(lambda) ~= [1 2]
   error('LAMBDA is not a scalar')
end

%  DATA

[nrows,ncols] = size(data);
if ncols > nrows
    tmp = ncols;
    ncols = nrows;
    nrows = tmp;
    data = data';
end
if ncols == 1
    data  = [data, (1:nrows)'];
    ncols = 2;
end
if ncols ~= 2
    error('DATA is neither a n-by-1 nor a n-by-2 array');
end

%  ---------------------------------------------------------------
%  Construct penalty matrix and 'b' vector for Ax=b.
%  ---------------------------------------------------------------

params   = getbasispar(basisobj);
numnodes = size(params.nodes,1);

nodeStruct.order     = params.order;
nodeStruct.nodes     = params.nodes;
nodeStruct.nodeindex = params.nodeindex;
nodeStruct.J         = params.J;
nodeStruct.metric    = params.metric;

%  ---------------------------------------------------------------
% construct mass matrix K0 
%  ---------------------------------------------------------------

K0 = mass(nodeStruct);

%  ---------------------------------------------------------------
% construct stiffness matrix K1
%  ---------------------------------------------------------------

K1 = stiff1(nodeStruct);

%  ---------------------------------------------------------------
% construct projection matrix on the space spanned by the columns of the 
% design matrix desmat
% ATTENZIONE: is it possible to get the projection matrix without 
% having to compute it as below here?
%  ---------------------------------------------------------------

if ~isempty(desmat)
    Q = qr(desmat, 0);
    H = Q * Q';
else
    H = [];
end

%  ---------------------------------------------------------------
% construct the block diagonal matrix L, having upper left block 
% given by I-H and zero otherwise
%  ---------------------------------------------------------------

inddata = data(:,1);
penalty = zeros(numnodes,1);
penalty(inddata) = 1;
indnodes = 1:numnodes;

%  ---------------------------------------------------------------
% construct vector b for system Ax = b
%  ---------------------------------------------------------------
    
b = sparse(zeros(numnodes*2,1));
b(inddata) = data(:,2);

L = sparse(indnodes,indnodes,penalty);
if ~isempty(desmat)
    L(inddata,inddata) = L(inddata,inddata) - H;
    b(1:numnodes,:)  = L * b(1:numnodes,:);
end

%  ---------------------------------------------------------------
% construct matrix A for system Ax=b.
%  ---------------------------------------------------------------
    
A  = [ L    -lambda*K1; ...
      K1            K0];

% solve system

bigsol = A\b;  

solution = bigsol(indnodes);
s = bigsol(indnodes+numnodes);
    
%  ---------------------------------------------------------------
% Make output objects
%  ---------------------------------------------------------------
    
smooth_fd  = fd(solution, basisobj);

if nargout > 1
    laplace_fd = fd(s, basisobj);
end