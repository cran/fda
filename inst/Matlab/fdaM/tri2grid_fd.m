function [uxy,tn,al2,al3] = tri2grid_fd(p,t,u,tn,al2,al3)
%TRI2GRID_FD Interpolate from triangular mesh to rectangular grid.
%
%       UXY = TRI2GRID(P,T,U,X,Y) computes the function values UXY
%       over the grid defined by the vectors X and Y, from the
%       function U with values on the triangular mesh defined by P and T.
%       Values are computed using linear interpolation in
%       the triangle containing the grid point. Vectors X and Y
%       must be increasing.
%
%       [UXY,TN,A2,A3] = TRI2GRID(P,T,U,X,Y) lists
%       the index TN of the triangle containing each grid point,
%       and interpolation coefficients A2 and A3.
%
%       UXY = TRI2GRID(P,T,U,TN,A2,A3), with TN, A2, and A3 computed
%       in an earlier call to TRI2GRID, interpolates using
%       the same grid as in the earlier call. This variant is,
%       however, much faster if several functions have to be
%       interpolated using the same grid.
%
%       For grid points outside of the triangular mesh, NaN is
%       returned in UXY, TN, A2, and A3.
%
%
%  This version is an adaption of function tri2grid in the pde toolbox
%  in Matlab by Jim Ramsay.  
%
%  Last modified 8 March 2011

small = 10000*eps;

if nargin==5,
    %  TRI2GRID(P,T,U,X,Y)  ... interpolation without known weights
  search = 1;
  x = tn;
  y = al2;
else
    %  TRI2GRID(P,T,U,TN,A2,A3) ... interpolation with known weights
  search = 0;
  ny = size(tn,1);
end

if search==1,

  nt = size(t,1);
  nx = length(x);
  ny = length(y);

  tt    = zeros(nt,3);
  tt(:) = p(t(:,1:3),1);
  xmin  = min(tt,[],2);
  xmax  = max(tt,[],2);
  tt(:) = p(t(:,1:3),2);
  ymin  = min(tt,[],2);
  ymax  = max(tt,[],2);

  %  x-indices of left vertices
  
  [xmin,ix] = sort(xmin);
  j = 1;
  for i = 1:nt,
    if j<=nx,
      while x(j)<xmin(i),
        j = j+1;
        if j>nx,
          break
        end
      end
    end
    xmin(i) = j;
  end
  xmin(ix) = xmin;

  %  x-indices of right vertices
  
  [xmax,ix] = sort(xmax);
  j = nx;
  for i = nt:-1:1,
    if j>=1,
      while x(j)>xmax(i),
        j = j-1;
        if j<1,
          break
        end
      end
    end
    xmax(i) = j;
  end
  xmax(ix) = xmax;

  %  y-indices of lower vertices
  
  [ymin,ix] = sort(ymin);
  j = 1;
  for i = 1:nt,
    if j<=ny,
      while y(j)<ymin(i),
        j = j+1;
        if j>ny,
          break
        end
      end
    end
    ymin(i) = j;
  end
  ymin(ix) = ymin;

  %  y-indices of upper vertices

  [ymax,ix] = sort(ymax);
  j = ny;
  for i = nt:-1:1,
    if j>=1,
      while y(j)>ymax(i),
        j = j-1;
        if j<1,
          break
        end
      end
    end
    ymax(i) = j;
  end
  ymax(ix) = ymax;

  %  set up coefficients for computing barycentric coordinates
  
  tricoef = tricoefCal(p, t);

  %  compute interpolation coefficients

  tn  = NaN*ones(ny,nx);
  al2 = tn;
  al3 = tn;
  %  run through triangles from left to right
  for i = 1:nt,
    if xmin(i)<=xmax(i) && ymin(i)<=ymax(i),
      a2 = p(t(i,2),:)-p(t(i,1),:);
      a3 = p(t(i,3),:)-p(t(i,1),:);
      b2 = [a3(2); -a3(1)];
      b3 = [a2(2); -a2(1)];
      b2 = b2./(a2*b2);
      b3 = b3./(a3*b3);
      %  run through mesh indices within triangle interval
      for j = xmin(i):xmax(i),
          for k = ymin(i):ymax(i),
              %  check to see if point (x(j),y(i)) is already
              %  assigned
              %  check if point (x(j),y(i)) is within triangle i
              ind = insideIndex(x(j), y(k), p, t(i,:), tricoef(i,:));
              %  if not already assigned and point is within
              %  triangle, interpolate data
%               disp([j,k,x(j),y(k),ind])
              if ~isnan(ind) && isnan(tn(k,j))
%               if ~isnan(ind)
                  r1p = [x(j);y(k)]-p(t(i,1),:)';
                  d2 = b2'*r1p;
                  if d2>=-small && d2<=1+small,
                      d3 = b3'*r1p;
                      if d3>=-small && d2+d3<=1+small,
                          tn(k,j)  = i;
                          al2(k,j) = d2;
                          al3(k,j) = d3;
                      end
                  end
                  % end of interpolation
              end
          end
      end
    end
        %  end of loop thourgh mesh indices
  end
  %  end of loop through triangles
end 

%  interpolation with given or computed coefficients

ind = ~isnan(tn);
uxy = NaN*ones(size(tn));
if ~isempty(ind)
    uxy(ind) = (1-al2(ind)-al3(ind)).*u(t(tn(ind),1)) + ...
                           al2(ind) .*u(t(tn(ind),2)) + ...
                           al3(ind) .*u(t(tn(ind),3));
end

