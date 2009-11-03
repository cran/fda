function integhat = sfdint(space_basis, v, integorder, nderivs)
%  Approximate the integral of a spatial functional data object over a 
%  triangular region
%
%  Arguments:
%  SPACE_BASIS ...  A spatial basis object of the FEM class
%  V           ...  A 3 by 2 matrix of x-y coordinates of vertices
%  INTEGORDER  ...  The order of approximation:  the approximation is exact
%                   for spatial polynomials of that order
%  NDERIVS     ...  A vector of length 2 containing the orders of derivatives
%                   for X and Y, respectively
%
%  Last modified 27 September 2011

if nargin < 4, nderivs = zeros(1,2);  end
if nargin < 3, integorder = 2;        end

%  specify the barycentric coordinates of quadrature points

switch integorder
    case 0
        B = ones(3,1)/3;
        W = 1;
    case 1
        B = eye(3);
        W = ones(3,1)/3;
    case 2
        B = [1, 1, 0; ...
             1, 0, 1; ...
             0, 1, 1]/2;
        W = ones(3,1)/3;
    case 3
        B = [[1, 0, 0; ...
              0, 1, 0; ...
              0, 0, 1]; ...
             [1, 1, 0; ...
              1, 0, 1; ...
              0, 1, 1]/2;
             [1, 1, 1]/3];
        W = [3; 3; 3; 8; 8; 8; 27]/60;
    case 4
        a = (3 + sqrt(3))/3;
        b = 2 - a;
        c = 2/3;
        B = [1, 1, 0; ...
             1, 0, 1; ...
             0, 1, 1; ...
             a, b, 0; ...
             b, a, 0; ...
             a, 0, b; ...
             b, 0, a; ...
             0, a, b; ...
             0, b, a; ...
             c, c, c]/2;
         W = [-1; -1; -1; 6; 6; 6; 6; 6; 6; 27]/60;
    otherwise
        error('Invalid value for ORDER.');
end
         
%  compute the quadrature points

xymat = B*v;         

%  evaluate the function at the quadrature points

evalvec = eval_FEM_fd(xymat(:,1),xymat(:,2), space_basis, nderivs);

%  approximate the integral

integhat = sparse(W'*evalvec);

             