function [newfdobj,laplacefd] = smooth_FEM_fd_Covar(data,desmat,fdobj,lambda)
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
% Last modified on 29 June 2010 by laura
% from original smooth_FEM_fd code last modified on 14 June 2010 bi Jim Ramsay.

%  DESMAT   the design matrix

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

% Projection matrix
% ATTENZIONE: is it possible to get the projection matrix without 
% having to compute it as below here?

H= desmat * (( desmat' * desmat )^(-1)) * desmat';

penalty = zeros(n,1);
penalty(data(:,1)) = 1;
bigH = sparse(1:n,1:n,penalty);
bigH(data(:,1),data(:,1))=bigH(data(:,1),data(:,1))-H;

b = sparse(zeros(n*2,1));
b(data(:,1)) = data(:,2);

bigHb  = sparse(zeros(n*2,1));
bigHb(1:n,:)  = bigH * b(1:n,:);


% construct matrix A for system Ax=b.

K0 = stiff0(nodemesh, nodes);
K2 = stiff2(nodemesh, nodes);
A  = [bigH -lambda*K2; K2 K0];

% solve system

bigsol = A\bigHb;
u = bigsol(1:n);
s = bigsol((n+1):(2*n));
newfdobj  = fd(u, basisobj);
laplacefd = fd(s, basisobj);
   