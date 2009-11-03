function K1 = stiff1(nodeStruct)

%STIFF1 produces the nnod*nnod stiffness matrix K1
%defined (K1)jk = int(dpsik/da*dpsij/da + dpsik/db*dpsij/db).
%
%Input: NODESTRUCT is a struct object produced by function makenodes.
%    It contains:
%        ORDER     ... The order of the element (1 or 2)
%        NODES     ... Coordinates of node points
%        NODEINDEX ... indices of node points for each element
%        JVEC      ... Jacobian of the affine transformation of each
%                      element to the master element
%        METRIC    ... The crossproduct of the inverse of the linear
%                      part of the transformation
%
%Output: K1 is an nnod*nnod matrix out which is
%        the sum of the nele element stiffness matrices
%        and the penalty stiffness matrix.
%        These i'th element matrix has (ij)'th element defined
%        as follows:  
%        Let psita and psitb be the partial derivatives of the
%        t'th shape function with respect to a and b (1<=t<=6).
%        Then the integral of the sum of products
%        (psija*psika+psijb+psikb) over the i'th element is
%        computed.  Then that value is assigned to the
%        (nodeindex(i,j),nodeindex(i,k))'th entry of the i'th elemental 
%        stiffness matrix and the other elements are given the value zero.
%
% Last modified on 9 May 2011 by Jim.

%  retrieve arrays from nodeStruct

order     = nodeStruct.order;
nodes     = nodeStruct.nodes;
nodeindex = nodeStruct.nodeindex;
Jvec      = nodeStruct.J;
metric    = nodeStruct.metric;

nele = size(nodeindex,1);
nnod = size(nodes,1);
K1   = sparse(nnod,nnod);

%  assemble the stiffness matrix

if order == 2
    
    %  values of K1 for master elements
    
    KXX = [ 3, -4,  1,  0,  0,  0; ...
           -4,  8, -4,  0,  0,  0; ...
            1, -4,  3,  0,  0,  0; ...
            0,  0,  0,  8,  0, -8; ...
            0,  0,  0,  0,  0,  0; ...
            0,  0,  0, -8,  0,  8]./6;

    KXY = [ 3,  0,  0,  0,  1, -4; ...
           -4,  4,  0, -4,  0,  4; ...
            1, -4,  0,  4, -1,  0; ...
            0, -4,  0,  4,  4, -4; ...
            0,  0,  0,  0,  0,  0; ...
            0,  4,  0, -4, -4,  4]./6;

    KYY = [ 3,  0,  0,  0,  1, -4; ...
            0,  8,  0, -8,  0,  0; ...
            0,  0,  0,  0,  0,  0; ...
            0, -8,  0,  8,  0,  0; ...
            1,  0,  0,  0,  3, -4; ...
           -4,  0,  0,  0, -4,  8]./6;

    %  assemble the stiffness matrix
    
    for el=1:nele
        ind    = nodeindex(el,:);
        K1M = (metric(el,1,1).*KXX  + metric(el,1,2).*KXY + ...
               metric(el,2,1).*KXY' + metric(el,2,2).*KYY);
        K1(ind,ind) = K1(ind,ind) + K1M.*Jvec(el);
    end;
    
elseif order == 1
    
    KXX = [ 1, -1,  0; ...
           -1,  1,  0; ...
            0,  0,  0]./2;

    KXY = [ 1,  0, -1; ...
           -1,  0,  1; ...
            0,  0,  0]./2;

    KYY = [ 1,  0, -1; ...
            0,  0,  0; ...
           -1,  0,  1]./2;
       
    %  assemble the stiffness matrix
    
    for el=1:nele
        ind      = nodeindex(el,:);
        K1M = (metric(el,1,1).*KXX  + metric(el,1,2).*KXY + ...
               metric(el,2,1).*KXY' + metric(el,2,2).*KYY);
        K1(ind,ind) = K1(ind,ind) + K1M.*Jvec(el);
    end;
    
else
    error('ORDER not 1 or 2.');
end




