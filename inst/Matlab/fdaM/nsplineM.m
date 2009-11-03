function [nsplinemat,P] = nsplineM(x, breaks, norder, nderiv, sparsewrd)
%  NSPLINEM  Computes values or derivative values of an N-spline basis
%  functions as well as a projection matrix to convert the B-spline
%  coefficients into N-spline coefficients (via linear projection) 
%  Arguments:
%  X         ... Argument values for which function values are computed
%  BREAKS    ... Increasing knot sequence spanning argument range
%  NORDER    ... Order of N-spline (one greater than degree) max = 19
%                Default 4.
%  NDERIV    ... Order of derivative required, default 0.
%  SPARSEWRD ... if 1, return in sparse form
%  Return:
%  NSPLINEMAT ... length(X) times number of basis functions matrix
%                 of Nspline values
%  P          ... projection matrix to convert B-spline coefficients to
%                 corresponding N-spline coefficients.

%  added by Kris Villez in August 2011 based on bsplinepen file in the
%  FDA toolbox by Jim Ramsay

%  last modified 31 October 2011 by Kris Villez: changed terminology

%  check dimensions of X and set up as a row vector
 
sizex = size(x);
if sizex(1) > 1 && sizex(2) > 1
    error('Argument X is not a vector.');
end
x = x(:);

%  set default argument values

if nargin < 5, sparsewrd = 1;  end
if nargin < 4, nderiv    = 0;  end
if nargin < 3, norder    = 4;  end

ndcon       =   2   ; % constraint at derivative 2

% first compute B-splines second derivative coefficients at first and last knot.
bsplinemat  =   bsplineM(breaks, breaks, norder, ndcon, sparsewrd)   ;

CoeffL  =   bsplinemat(1,bsplinemat(1,:)~=0)                ;   % left side constraint coeff
CoeffR  =   rot90(bsplinemat(end,bsplinemat(end,:)~=0),2)   ;   % right side constraint coeff
nL      =   length(CoeffL)          ;
nR      =   length(CoeffR)          ;

% projection matrix
nb      =   length(breaks) + norder - 2     ;

PL      =   sparse(eye(nL))         ;
PL(1,:) =   -CoeffL(:)/CoeffL(1)    ;
PL      =   PL(:,2:end)             ;
PR      =   sparse(eye(nR))         ;
PR(1,:) =  -CoeffR(:)/CoeffR(1)     ;
PR      =   PR(:,2:end)             ;

P           =   zeros(nb,nb-2)  ;
vr          =   1:nL            ;
vc          =   1:nL-1          ;
P(vr,vc)    =   PL              ;
vr          =   nL+1:nb-nR      ;
vc          =   nL:nb-nR-1      ;
P(vr,vc)    =   eye(length(vc)) ;
vr          =   nb+1-(1:nR)     ;
vc          =   nb-1-(1:nR-1)   ;
P(vr,vc)    =   PR              ;

% now compute B-splines desired values/derivative.
bsplinemat  =   bsplineM(x, breaks, norder, nderiv, sparsewrd)   ;

nsplinemat  =   bsplinemat*P    ;

if sparsewrd, nsplinemat = sparse(nsplinemat); end
