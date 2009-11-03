function [MS_amp, MS_pha, RSQR, C] = AmpPhaseDecomp(xfd, yfd, hfd, rng)
%  Computes the amplitude-phase decomposition for a registration.

%  Arguments:
%  XFD  ...  FD object for unregistered functions
%  YFD  ...  FD object for registered functions
%  Hfd  ...  FD object for warping functions
%  RNG  ...  Sub-interval over which the decomposition is computed

%  Returns:
%  MS_amp ... mean square for amplitude variation 
%  MS_pha ... mean square for amplitude variation 
%  RSQR   ... squared correlation measure of prop. phase variation 
%  C      ... constant C

%  Last modified 21 April 2009

%  extract information from XFD

xbasis  = getbasis(xfd);
xrng    = getbasisrange(xbasis);

%  set default interval RNG

if nargin < 4, rng = xrng;  end

%  check RNG

if rng(1) < xrng(1) || rng(2) > xrng(2)
    error('RNG values outide XRNG interval.');
end
if rng(1) >= rng(2)
    error('RNG values not strictly increasing.');
end

%  set up a fine mesh of values for numerical integration

nxbasis = getnbasis(xbasis);
nfine   = max([201,10*nxbasis]);
tfine   = linspace(rng(1),rng(2),nfine)';
delta   = tfine(2) - tfine(1);

%  evaluate arguments at fine mesh values

xfine   = eval_fd(tfine, xfd);
yfine   = eval_fd(tfine, yfd);
Dhfine  = eval_fd(tfine, hfd, 1);

%  means of unregistered and registered functions

mufine  = mean(xfine,2);
etafine = mean(yfine,2);

%  integrated squared mean values

intetasqr = delta.*trapz(etafine.^2);
intmusqr  = delta.*trapz( mufine.^2);

%  covariance between Dh and y values

covDhSy = zeros(nfine,1);
for i=1:nfine
    Dhi        = Dhfine(i,:);
    Syi        = yfine(i,:).^2;
    Covmat     = cov(Dhi', Syi');    
    covDhSy(i) = Covmat(1,2);
end

%  integrated covariance value

intcovDhSy = delta.*trapz(covDhSy);

%  integrated squared registered and registered residual values

N = size(xfine,2);
intysqr = zeros(N,1);
intrsqr = zeros(N,1);
rfine   = yfine - etafine*ones(1,N);
for i=1:N
    intysqr(i) = delta.*trapz(yfine(:,i).^2);
    intrsqr(i) = delta.*trapz(rfine(:,i).^2);
end

%  compute results

C      = 1 + intcovDhSy/mean(intysqr);
MS_amp = C*mean(intrsqr);
MS_pha = C*intetasqr - intmusqr;
RSQR   = MS_pha/(MS_amp+MS_pha);
