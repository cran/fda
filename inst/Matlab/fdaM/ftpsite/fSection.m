function fdout = fSection(fdin, cutpoints, basisobj, n, ...
                          lambda, Lfdin, Lfdout)
%  fSection sections an fdobject FDIN by representing the
%  input function restricted to each of a set of
%  intervals as a separate fd object using the basis
%  supplied in argument BASISOBJ.  
%  This procedure implies that each of the section objects
%  is defined over a fixed common range that is
%  supplied in the basis object.  This common range need
%  not correspond in any way to the range of the original
%  fdobject, and could, for example, be simply [0,1].
%  If there are M sections, each is defined by M-1 strictly
%  increasing cutpoints in argument CUTPOINTS. The first
%  of these cutpoints must be greater than the lower limit
%  in the range over which FDIN is defined, and the last
%  must be strictly less that the upper limit for FDIN. 
%  FDOUT will contain M replications.
%  Each section's fdobj will be defined by smoothing
%  N equally spaced values over the section, where 
%  N >= the number of basis functions in BASISOBJ.
%  N may be supplied as an argument, or it may be default
%  to the number of basis functions in BASISOBJ.
%  The level of smoothing is controlled by smoothing
%  parameter lambda, which defaults to 1e-10; and
%  by linear differential operator LFDIN, which
%  defaults to 2.  If a warning message about singularity
%  is issued by smooth_basis, lambda may need to be 
%  increased, or, alternatively, n can be increased.  
%  This interpolation procedure will guarantee that the
%  position and the shape of a number of derivatives will
%  be preserved over the interval, where the number of
%  derivatives to be preserved depends on the order of 
%  BASISOBJ.  
%  The sections can correspond to one or more derivatives of
%  FDIN by controlling the argument LFDOUT, which defaults
%  to 0.  Of course, the order of the highest derivative
%  in LFDOUT must be consistent with the order of BASISOBJ.  

%  Last modified 5 April 2007

%  set default value for NDERIV

if nargin < 8, Lfdout = int2Lfd(0);  end
if nargin < 7, Lfdin  = int2Lfd(2);  end
if nargin < 6, lambda = 1e-10;       end

%  Get basis information from FDIN

basisin = getbasis(fdin);
rangein = getbasisrange(basisin);

%  get basis information for sections

nbasis = getnbasis(basisobj);
rangeout = getbasisrange(basisobj);

%  set default value for N

if nargin < 4, n = nbasis;  end
argvals = linspace(rangeout(1),rangeout(2),n)';

%  Check the cutpoints

M = length(cutpoints) + 1;
if any(diff(cutpoints) <= 0)
    error('Cutpoints are not strictly increasing.');
end
if min(cutpoints) <= rangein(1) || ...
   max(cutpoints) >= rangein(2)
    error('One or more cutpoints are out of range.');
end

%  Add lower and upper limits to cutpoints

cutpoints = cutpoints(:);  %  ensures column vector
breaks    = [rangein(1); cutpoints; rangein(2)];

%  loop through intervals

ymat = zeros(nbasis,M);
for i=1:M
    rangei = [breaks(i), breaks(i+1)];
    xi = linspace(rangei(1),rangei(2),n)';
    yi = eval_fd(xi, fdin, Lfdout);
    ymat(:,i) = yi;
end

fdParobj = fdPar(basisobj, Lfdin, lambda);
fdout = smooth_basis(argvals, ymat, fdParobj);



