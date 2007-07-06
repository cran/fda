function plot(bifdParobj, nx)
%  plots a bivariate functional parameter object
%  Arguments:
%  BIFDOBJ ... A bivariate functional data object to be plotted.
%  MATPLOT ... If MATPLT is nonzero, all curves are plotted in a 
%              single plot.
%              Otherwise, each curve is plotted separately, and the
%              next curve is plotted when the mouse is clicked.
%  HREF    ... If HREF is nonzero, a horizontal dotted line through 0 
%              is plotted.
%  NX      ... Number of equally spaced argument values to use for 
%              plotting.  Default is 101.  

%  Last modified 14 October 2003

%  set default arguments
if nargin < 2, nx = 101;   end

%  plot the functional data object

plot(bifdParobj.fd, nx)

%  plot the object to which Lfd has been applied

display('Press any key to see LFDOBJ applied to the parameter.);

plot(bifdParobj.fd, bifdParobj.Lfd, matplt, href, nx)
