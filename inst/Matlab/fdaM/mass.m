function K0 = mass(nodeStruct)

%MASS_NEW produces the mass matrix containing integrals of products of
%  nodal functions.  
%
%Input: NODESTRUCT is a struct object produced by function makenodes_new.
%    It contains:
%        ORDER     ... The order of the element (1 or 2)
%        NODES     ... Coordinates of node points
%        NODEINDEX ... indices of node points for each element
%        JVEC      ... Jacobian of the affine transformation of each
%                      element to the master element
%
%Output: K0: the NNOD by NNOD matrix of sums of products of nodal basis
%        functions.
%        For each element i, the integral of the product 
%        of the j'th and k'th shape functions over the i'th element is
%        computed.  Then that value is the 
%        (NODEINDEX(i,j),NODEINDEX(i,k))'th entry of the i'th elemental 
%        mass matrix.
%
% Last modified on 27 September 2011 by Jim.

%  retrieve arrays from nodeStruct

order     = nodeStruct.order;
nodes     = nodeStruct.nodes;
nodeindex = nodeStruct.nodeindex;
Jvec      = nodeStruct.J;

nele  = size(nodeindex,1);
nnod  = size(nodes,1);

if order ==2
    
    %  the integrals of products of basis functions for master element:
    
    K0M = [[ 6  0 -1 -4 -1  0]; ...
           [ 0 32  0 16 -4 16]; ...
           [-1  0  6  0 -1 -4]; ...
           [-4 16  0 32  0 16]; ...
           [-1 -4 -1  0  6  0]; ...
           [ 0 16 -4 16  0 32]]./360;
    
    %  assemble the mass matrix
    
    K0 = sparse(nnod,nnod);
    for el=1:nele
        ind = nodeindex(el,:);
        K0(ind,ind) = K0(ind,ind) + K0M.*Jvec(el);
    end;
    
elseif order == 1
    
    %  the integrals of products of basis functions for master element:
    
    K0M = [[ 2  1  1]; ...
           [ 1  2  1]; ...
           [ 1  1  2]]./24;
    
    %  assemble the mass matrix
    
    K0 = sparse(nnod,nnod);
    for el=1:nele
        ind = nodeindex(el,:);
        K0(ind,ind) = K0(ind,ind) + K0M.*Jvec(el);
    end;
    
else
    error('ORDER not 1 or 2.');
end
