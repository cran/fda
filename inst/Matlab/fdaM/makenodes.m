function [nodes, nodemesh] = makenodes(p,t)

%MAKENODES produces a matrix NODES of coordinates for all of the
%nodes to be used and a matrix NODEMESH defining which nodes
%correspond to each element.  The midpoint of each edge is
%computed and added to P and the row index of that midpoint
%is then added to the rows of T containing that edge.

%Input: P is an NVER by 2 matrix containing the X and Y
%	coordinates of each of the nvert p in the
%	right-hand rule mesh.
%   The call can use P directly (see below).
%	T contains the indices of the points in P defining the triangles in
%   counter-clockwise order.  T can have either 3 or 4 columns.  The
%   fourth column is ignored.
%
%Output: NODES:  a NNODES*2 matrix whose i'th row contains
%		the coordinates of the i'th nodal variable.
%       Nodes for the second order Lagrange element consist of vertices 
%       and midpoints of edges, that is, 6 per triangle.
%
%		The first NVER rows of NODES is P, and the remainder are
%       the edge midpoints.
%
%	 NODEMESH:  an NNODES*6 matrix whose i'th row
%		contains the row numbers (in NODES) of the
%		nodal variables defining the i'th finite 
%		element.  If the i'th row of T is [V1 V2 V3]
%		then the i'th row of NODEMESH is
%		[V1 V(12) V2 V(23) V3 V(31)], where Vi is the
%		row number of the i'th point and V(ij) is the 
%		row number of the midpoint of the edge defined
%		by the i'th and j'th points.
%
% Last modified on June 14, 2010 by Jim Ramsay

%  The first rows of nodes are the vertices

nodes = p;

nele = size(t,1);
nver = size(p,1);

ind  = [ 1 2 ; 2 3 ; 3 1 ];

if nargout > 1
    nodemesh = sparse(nele,6);
    nodemesh(:,[1 3 5]) = t(:,1:3);
    rec  = sparse(nver,nver);
end

for i = 1:nele
  for j = 1:3
    if rec(t(i,ind(j,1)),t(i,ind(j,2)))==0
      nodes = [nodes; [.5 .5]*nodes(t(i,ind(j,:)),:)];
      if nargout > 1
          rec(t(i,ind(j,1)),t(i,ind(j,2))) = size(nodes,1);
          rec(t(i,ind(j,2)),t(i,ind(j,1))) = size(nodes,1);
          nodemesh(i,2*j) = size(nodes,1);
      end
    else
        if nargout > 1
            nodemesh(i,2*j) = rec(t(i,ind(j,1)),t(i,ind(j,2)));
        end
    end;
  end;
end;

