function [felsplobj,laplacefd] = smooth_FEM_fd_new(data,fdobj,lambda)
% SMOOTH_FELSPLINE Compute a new solution for a FELspline
%  problem 
%
%     Arguments:
% FELSPLOBJ a FELspline object, constructed by SOLVE_FELSPLINE.
% LAMBDA    a scalar smoothing parameter
% DATA      (optional) a n-by-2 new set of observations.  DATA(:,1)
%           indexes the points (FELSPLOBJ.POINTS) at which the 
%           values in DATA(:,2) were observed.
%
%     Output:
% FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
%                 vector resulting from smoothing
% LAPLACEFD  ...  A FD object of the FEM type for the value of the 
%                 Laplace operator if order == 2, or empty if order == 1
%
% Last modified on 26 August 2010.

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

basisobj  = getbasis(fdobj);
params    = getbasispar(basisobj);
p         = params.p;
e         = params.e;
t         = params.t;
nodes     = params.nodes;
nodeindex = params.nodeindex;
numnodes  = size(nodes,1);
order     = params.order;

nodeStruct.order     = order;
nodeStruct.nodes     = nodes;
nodeStruct.nodeindex = params.nodeindex;
nodeStruct.J         = params.J;
nodeStruct.metric    = params.metric;

%  ---------------------------------------------------------------
% construct mass matrix K0 
%  ---------------------------------------------------------------

K0 = mass_new(nodeStruct);

%  ---------------------------------------------------------------
% construct stiffness matrix K0 and K1.
%  ---------------------------------------------------------------

K1 = stiff1_new(nodeStruct);

%  ---------------------------------------------------------------
%              second order elements
% construct the penalty matrix P with ones on diagonal at data points
%  ---------------------------------------------------------------

penalty = zeros(numnodes,1);
penalty(data(:,1)) = 1;
indnodes = 1:numnodes;
P = sparse(indnodes,indnodes,penalty);

%  ---------------------------------------------------------------
% construct the solution
%  ---------------------------------------------------------------

if order == 2
    
    %  ---------------------------------------------------------------
    % construct vector f for system Ax=f
    %  ---------------------------------------------------------------
    
    f = sparse(zeros(numnodes*2,1));
    f(data(:,1)) = data(:,2);
    
    %  ---------------------------------------------------------------
    % construct matrix A for system Ax=f.
    %  ---------------------------------------------------------------
    
    A = [P  -lambda*K1; ...
        K1         K0];
    
    % solve system
    
    bigsol = A\f;
    
    solution = bigsol(indnodes);
    s = bigsol(indnodes+numnodes);
    
    % Make FELspline object
    
    felsplobj = fd(solution, basisobj);
    laplacefd = fd(s, basisobj);
    
elseif     order == 1
    
    %  ---------------------------------------------------------------
    %              first order elements
    % construct the penalty matrix P with ones on diagonal at data points
    %  ---------------------------------------------------------------
    
    f = zeros(numnodes,1);
    f(data(:,1)) = data(:,2);
    
    solution = (P + lambda.*K1)\f;
    
    felsplobj = fd(solution, basisobj);
    laplacefd = [];
    
else
    error('NDERIV is neither 1 nor 2.');
end

%  ------------------------------------------------------------------------

function fm = rhs(p,t,f,nderiv)
%RHS Assembles right hand side in a PDE problem.

%  modified 8 June 2010 by Jim

%  set default value for a

if nargin < 3, nderiv = 1;  end

if nderiv == 1
    
np = size(p,2); % Number of points

% Corner point indices

it1 = t(1,:);
it2 = t(2,:);
it3 = t(3,:);

% Triangle geometries:

ar = pdetrg(p,t);

% Find midpoints of triangles

x = (p(1,it1)+p(1,it2)+p(1,it3))/3;
y = (p(2,it1)+p(2,it2)+p(2,it3))/3;

% ---------------------------  RHS  ---------------------------------------

  nrf = size(f,1);

  if nrf==1, % Scalar f
    fm = pdeasmf(it1,it2,it3,np,ar,x,y,f);
  else
    error('Wrong number of rows of f.');
  end % size(f,1)

elseif nderiv==2
    
    fm = f(:);
    
else
        error('NDERIV is neither 1 nor 2.');
end

   