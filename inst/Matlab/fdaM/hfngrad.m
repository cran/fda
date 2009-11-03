function [hfine, dhdc, dhdt, d2hdtdc] = ...
                     hfngrad(xfine, xlo, xhi, Wfd, basiscell)
% HFNGRAD computes the values at times in XFINE of:
%    the warping h
%    the gradient of h with respect to coefficients
%  h(t) is defined as h(t) = M(t) / M(XHI-XLO) where
%                     M(t) = int_0^t exp[W(u)] du
%  It is assumed that the last value in XFINE = XHI - XLO = WIDTH.

%  Last modified 18 April 2011 by Jim Ramsay

nfine    = length(xfine);
width    = xhi - xlo;
%  evaluate function M at points in XFINE
Mfine    =   monfn(xfine, Wfd, basiscell);
%  evaluate the gradient of M at points in XFINE
dMdcfine = mongrad(xfine, Wfd, basiscell);
%  get normalizing constant = max(Mfine).
Mmax     = Mfine(nfine);
%  get normalizing gradient elements
dMdcmax  = dMdcfine(nfine,:);
%  define HFINE
hfine    = xlo + width.*Mfine./Mmax;
%  define DHDC
dhdc     = width.*(Mmax.*dMdcfine - Mfine*dMdcmax)./Mmax^2;
%  ensure that minimum and maximum elements of HFINE are correct
hfine(1)     = xlo;
hfine(nfine) = xhi;
if nargout > 2
    %  evaluate dhdt
    dMdt = exp(eval_fd(xfine, Wfd));
    dhdt = width.*dMdt/Mmax;
    %  evaluate gradient of dhdt
    Wbasis  = getbasis(Wfd);
    nbasis  = getnbasis(Wbasis);
    d2Mdtdc = eval_basis(xfine, Wbasis).*repmat(dMdt,1,nbasis);
    d2hdtdc = width.*(Mmax.*d2Mdtdc - dMdt*dMdcmax)./Mmax^2;
else
    dhdt    = [];
    d2hdtdc = [];
end
