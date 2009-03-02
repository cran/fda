function yshift = shifty(x, y, shift)
%SHIFTY estimates value of Y for periodic data for
%       X shifted by amount SHIFT.
%  It is assumed that X spans interval over which functionis periodic.
%  Last modified 30 December 2000

ydim = size(y);
nrep = ydim(2);
if length(ydim) <= 2
   nvar = 1;
else
   nvar = ydim(3);
end

if shift == 0
   yshift = y;
   return;
end

n   = ydim(1);
xlo = min(x);
xhi = max(x);
wid = xhi - xlo;
if shift > 0
   while shift > xhi,  shift = shift - wid;  end
   ind = 2:n;
   x2 = [x; x(ind)+wid];
   xshift = x + shift;
   if nvar == 1
      y2 = [y; y(ind,:)];
      yshift = interp1(x2, y2, xshift);
   else
      yshift = zeros(n,nrep,nvar);
      for ivar=1:nvar
         y2 = [squeeze(y(:,:,ivar)); squeeze(y(ind,:,ivar))];
         yshift(:,:,ivar) = interp1(x2, y2, xshift);
      end
   end
else
   while shift < xlo - wid, shift = shift + wid;  end
   ind = 1:n-1;
   x2 = [x(ind)-wid; x];
   xshift = x + shift;
   if nvar == 1
      y2 = [y(ind,:);   y];
      yshift = interp1(x2, y2, xshift);
   else
      yshift = zeros(n,nrep,nvar);
      for ivar=1:nvar
         y2 = [squeeze(y(ind,:,ivar)); squeeze(y(:,:,ivar))];
         yshift(:,:,ivar) = interp1(x2, y2, xshift);
      end
   end
end
