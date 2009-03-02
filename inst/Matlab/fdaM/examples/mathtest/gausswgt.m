function [wgt, integvec] = gausswgt(thetaq, nbasis, norder)
%  Quadrature weights for B-spline test functions

%  set up fine mesh for estimating integrals

nq        = length(thetaq);
thetarng  = [min(thetaq), max(thetaq)];
delta     = 1e-3;
thetafine = min(thetaq):delta:max(thetaq);
nfine     = length(thetafine);

%  set up spline basis

norder   = 4;
nbasis   = nq + norder - 2;
basisobj = create_bspline_basis(thetarng, nbasis);

basisnq  = getbasismatrix(thetaq,    basisobj);
basisfin = getbasismatrix(thetafine, basisobj);

%  compute std. normal density values

kernelfn = (exp(-thetafine.^2/2)/sqrt(2*pi))';

%  estimate integrals by trapezoidal rule

integvec = delta.*(basisfin'*kernelfn - 0.5.*(basisfin(1,:)'.*kernelfn(1) + basisfin(nfine,:)'.*kernelfn(nfine)));

%  solve for weights giving least squares solution

wgt = basisnq'\integvec;
