function newfdobj = smooth_FEM_simple(data,fdobj,lambda)
% SMOOTH_FELSPLINE Compute a new solution for a FELspline
%  problem for which the stiffness matrices have already 
%  been constructed.  Note, if stiffness matrices have not 
%  yet been constructed, use SOLVE_FELSPLINE.
% 1)SMOOTH_FELSPLINE(FELSPLOBJ,LAMBDA) smooths the data already
%      contained in the FELspline object, FELSPLOBJ.
% 2)SMOOTH_FELSPLINE(FELSPLOBJ,LAMBDA,DATA) smooths new data on 
%      on the triangulation used to build the FELspline object.
%
%     Arguments:
% FELSPLOBJ a FELspline object, constructed by SOLVE_FELSPLINE.
% LAMBDA    a scalar smoothing parameter
% DATA      (optional) a n-by-2 new set of observations.  DATA(:,1)
%           indexes the points (FELSPLOBJ.POINTS) at which the 
%           values in DATA(:,2) were observed.
%
%     Output:
% NEWFDOBJ  ...  A FD object of the FEM type defined by the coefficient
%   vector resulting from smoothing
%
% Last modified on 25 June 2010.

%  check arguments

if ~isa_fd(fdobj)
   error('FDOBJ is not a FD object');
end

if  ~isa(lambda,'double')
   error('LAMBDA is not numeric')
elseif size(lambda) ~= [1 2]
   error('LAMBDA is not a scalar')
end

%  check data argument

if nargin<3
   data=getdata(fdobj);
elseif size(data,2)~=2
   if size(data,1)~=2
      error('DATA is not a n-by-2 array')
   else
      data=data';
   end
end

%  Construct penalty matrix and 'b' vector for Ax=b.

basisobj = getbasis(fdobj);
params   = getbasispar(basisobj);
nodes    = params.nodes;
nodemesh = params.nodemesh;

n = size(nodes,1);
penalty = zeros(n,1);
penalty(data(:,1)) = 1;
P = sparse(1:n,1:n,penalty);

% construct matrix A for system Ax=b.

K2 = stiffGrad(nodemesh, nodes);

A = P + lambda*K2;
b = zeros(n,1);
b(data(:,1)) = data(:,2);

% solve system

bigsol = A\b;
u = bigsol(1:n);

newfdobj  = fd(u, basisobj);
   