function nodeStruct = makenodes(p, t, order)

%  MAKENODES computes information about a triangulation of a polyhedral
%  domain and the polynomial basis functions associated with each
%  triangle.  A triangle combined with the specification of these 
%  bivariate basis polynomials is called a finite element.
%
%  Argument P contains the NVER vertices of the triangles in the
%  triangulation, and argument T specifies the three points in P defining
%  each of the NELE triangleS.  The degree of the bivariate basis functions  
%  is defined by argument ORDER. 
%
%  MAKENODES produces a struct object whose fields contain information 
%  about the points input in argument P and the triangles input in 
%  argument T conditional on the order of the Lagrangian element input
%  in argument ORDER.
%
%  Arguments:
%  P     ...  An NVER by 2 matrix of X-Y coordinates of each vertex
%             in the triangulation.  If P is input as a 2 by NVER
%             matrix, it is transposed
%  T     ...  An NELE by 3 matrix, each row of which contains the
%             indices of the points in P that define a triangle.  These
%             indices must be specified in counter-clockwise order
%             beginning with any of the three vertices.  If T is
%             input as a 3 by NELE matrix, it is transposed, and if
%             input as a 4 by NELE matrix, it's first three rows are
%             transposed to define T.
%  ORDER ...  The degree of the bivariate polynomials defining the
%             basis functions.  
%             If ORDER = 1, there are three basis functions.  
%             Each of these is  linear, equal to 1 on at one of the
%             vertices, and equal to zero on the edge of the triangle
%             opposite to this vertex.
%             If ORDER = 2, there are six basis functions.
%             each of these is quadratic, equal to 1 on either a vertex
%             or at a mid-point of an edge, and vanishes on edges not
%             containing the point at which it is 1.
%
%  Output:    A struct object with fields:
%  ORDER ...  the argument ORDER, which is the degree of each 
%             bivariate basis polynomial
%  NODES ...  a NNODES*2 matrix whose i'th row contains
%		      the coordinates of the i'th nodal variable.
%             NODES for the second order element consist of vertices and 
%             midpoints of edges, that is, 6 per triangle.
%		      The first NVER rows of NODES is P, and the remainder are
%             the edge midpoints.
%             Nodes for the first order element consist of only vertices 
%             in P.
%  NODEINDEX ... for ORDER == 2, an nele*6 matrix whose i'th row
%		      contains the row numbers (in NODES) of the
%		      nodal variables defining the i'th finite 
%		      element.  If the i'th row of P is [V1 V2 V3]
%		      then the i'th row of nodeindex is
%		      [V1 V(12) V2 V(23) V3 V(31)], where Vi is the
%		      row number in NODES of the i'th vertex and V(ij) is the 
%		      row number in NODES of the midpoint of the edge defined
%		      by the i'th and j'th points.
%             If NORDER == 1, NODEINDEX is T.
%  J     ...  Each triangle defines a linear transformation A from a fixed
%             right master traingle with vertices (0,0), (1,0) and (0,1)
%             to itself.  J contains one-half the Jacobian of this 
%             transformation, which is also the area of the triangle.
%  METRIC ... The 2 by 2 metric matrix inv(A)' inv(A)
%
% Last modified on 9 May 2011 by Jim Ramsay.

%  Set default order

if nargin < 3, order = 2;  end

if size(p,2)>2
    %  transpose if points and triangles come from pde toolbox
    nodes = p';
    t = t(1:3,:)';
else
    nodes = p;
end;

nele = size(t,1);  %  number of elements
nver = size(p,1);     %  number of vertices

Jvec   = zeros(nele,1);    %  vector of jacobian values
metric = zeros(nele,2,2);  %  3-d array of metric matrices

if order == 2
    
    rec = sparse(nver,nver);
    ind = [ 1 2 ; 2 3 ; 3 1 ];
    nodeindex = zeros(nele,6);
    nodeindex(:,[1 3 5]) = t;
    
    for i = 1:nele
        for j = 1:3
            if rec(t(i,ind(j,1)),t(i,ind(j,2)))==0
                nodes = [nodes; [.5 .5]*nodes(t(i,ind(j,:)),:)];
                rec(t(i,ind(j,1)),t(i,ind(j,2))) = ...
                    size(nodes,1);
                rec(t(i,ind(j,2)),t(i,ind(j,1))) = ...
                    size(nodes,1);
                nodeindex(i,2*j) = size(nodes,1);
            else
                nodeindex(i,2*j) = ...
                    rec(t(i,ind(j,1)),t(i,ind(j,2)));
            end;
        end;
        
        %  deviations of vertices 2 and 3 from vertex 1.  These
        %  define the 2 x 2 transformation matrix
        %  A = [diff1x  diff1y;  diff2x  diff2y]
        
        diff1x = nodes(nodeindex(i,3),1) - nodes(nodeindex(i,1),1);
        diff1y = nodes(nodeindex(i,3),2) - nodes(nodeindex(i,1),2);
        diff2x = nodes(nodeindex(i,5),1) - nodes(nodeindex(i,1),1);
        diff2y = nodes(nodeindex(i,5),2) - nodes(nodeindex(i,1),2);
        
        %  Half the Jacobian of the transformation from the
        %  master triangle to this triangle, which is also the area of 
        %  this triangle.
        
        Jvec(i) = (diff1x.*diff2y - diff2x.*diff1y)/2;
        
        %  Compute contravariant transformation matrix, also called
        %  the metric matrix.
        
        %  Compute the inverse transformation:
        
        Ael = [ diff2y -diff1y; ...
               -diff2x  diff1x]./Jvec(i);

        %  Compute metric matrix
        
        metric(i,:,:) = Ael'*Ael;
    end;
    
elseif order == 1
    
    nodeindex = t(:,1:3);
    
    for i=1:nele
        
        %  deviations of vertices 2 and 3 from vertex 1
        
        diff1x = nodes(nodeindex(i,2),1)-nodes(nodeindex(i,1),1);
        diff1y = nodes(nodeindex(i,2),2)-nodes(nodeindex(i,1),2);
        diff2y = nodes(nodeindex(i,3),2)-nodes(nodeindex(i,1),2);
        diff2x = nodes(nodeindex(i,3),1)-nodes(nodeindex(i,1),1);
        
        Jvec(i) = (diff1x.*diff2y - diff2x.*diff1y)/2;
        
        Ael = [ diff2y -diff1y; ...
               -diff2x  diff1x]./Jvec(i);

        metric(i,:,:) = Ael'*Ael;
    end
else
    error('ORDER not 1 or 2.');
end

nodeStruct.order     = order;
nodeStruct.nodes     = nodes;
nodeStruct.nodeindex = nodeindex;
nodeStruct.J         = Jvec;
nodeStruct.metric    = metric;