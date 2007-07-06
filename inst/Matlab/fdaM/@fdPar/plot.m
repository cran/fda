function plot(fdParobj, matplt, href, nx)
%  plots a functional parameter object
%  Arguments:
%  FDOBJ   ... A functional data object to be plotted.
%  MATPLOT ... If MATPLT is nonzero, all curves are plotted in a 
%              single plot.
%              Otherwise, each curve is plotted separately, and the
%              next curve is plotted when the mouse is clicked.
%  HREF    ... If HREF is nonzero, a horizontal dotted line through 0 
%              is plotted.
%  NX      ... Number of equally spaced argument values to use for 
%              plotting.  Default is 101.  

%  Last modified 9 September 2003

%  set default arguments
if nargin < 4, nx = 101;   end
if nargin < 3, href = 1;   end
if nargin < 2, matplt = 1; end

%  plot the functional data object

fdobj = fdParobj.fd;
plot(fdobj, int2Lfd(0), matplt, href, nx)

type = getbasistype(getbasis(fdobj));

if ~strcmp(type, 'const')
    
    %  plot the object to which Lfd has been applied
    
    display('Press any key to see LFDOBJ applied to the parameter.');
    
    pause;
    
    plot(fdParobj.fd, fdParobj.Lfd, matplt, href, nx)
    
end
