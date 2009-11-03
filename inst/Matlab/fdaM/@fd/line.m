function Hline = line(fdobj, Lfdobj)
%  LINE adds the curves in FD to an already existing plot.
%  Arguments:
%  FDOBJ  ... A functional data object
%  LFDOBJ  ... A linear differential operator object

%  last modified 3 March 2009

if nargin < 2, Lfdobj = int2Lfd(0);  end

if ~isa_fd(fdobj)
    error('FD is not a functional data object.');
end

nbasis = getnbasis(getbasis(fdobj));
nfine  = max([101,10*nbasis+1]); 

Lfdobj = int2Lfd(Lfdobj);
if ~isa_Lfd(Lfdobj)
    error ('LFD is not a linear differential operator object.');
end

coefd = size(getcoef(fdobj));
ndim  = length(coefd);
if ndim < 3
   rangex = getbasisrange(getbasis(fdobj));
   tfine  = linspace(rangex(1),rangex(2),nfine)';
   fdmat  = eval_fd(tfine, fdobj, Lfdobj);
   if nargout > 0
       Hline = line(tfine, fdmat);
   else
       line(tfine, fdmat);
   end
else
   error('Function LINE cannot be used with three dimensions.');
end

