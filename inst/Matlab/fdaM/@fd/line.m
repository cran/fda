function Hline = line(fd, Lfdobj, nx)
%  LINE adds the curves in FD to an already existing plot.
%  Arguments:
%  FD      ... A functional data object
%  LFDOBJ  ... A linear differential operator object
%  NX      ... The number of discrete values to be plotted.

%  last modified 27 November 2003

if ~isa_fd(fd)
    error('FD is not a functional data object.');
end

nbasis = getnbasis(getbasis(fd));

if nargin < 3, nx     = 10*nbasis+1; end
if nargin < 2, Lfdobj = int2Lfd(0);  end

Lfdobj = int2Lfd(Lfdobj);
if ~isa_Lfd(Lfdobj)
    error ('LFD is not a linear differential operator object.');
end

coefd = size(getcoef(fd));
ndim  = length(coefd);
if ndim < 3
   rangex = getbasisrange(getbasis(fd));
   x      = linspace(rangex(1),rangex(2),nx);
   fdmat  = eval_fd(x, fd, Lfdobj);
   Hline = line(x, fdmat);
else
   error('Too many dimensions');
end

