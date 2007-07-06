function [GCV, DGCV, df] = llbda2gcv(loglambda, Y, Z, R, factor)
% LLBDA2GCV computes the GCV criterion and its derivative
%  as a function of log_{10} \lambda where \lambda is the
%  smoothing parameter in a call to SMOOTH_BASIS
%  Arguments:
%  LOGLAMBDA ... logarithm to base 10 of \lambda
%  Y         ... data that are smoothed, n by N matrix where
%                n is the number of values per record that are
%                smoothed and N is the number of records.
%  Z         ... n by K matrix of basis functions where K
%                is the number of basis functions
%  R         ... K by K penalty matrix
%  FACTOR    ... Chong Gu's correction factor, suggested as 
%                1.2 or 1.4.  

%  Last modified 2 December 2006

%  In the following computation, no attempt is made to use
%  matrix decompositions since normally both Z and R will be
%  band-structured.

if nargin < 5, factor = 1;  end
if nargin < 4
    error('Less than four arguments.');
end

logtrR = log10(trace(R));

if logtrR + loglambda > 10
    warning('Wid:condition', ...
        'Condition number is too large for stable results.');
end

lambda = 10^loglambda;   %  lambda

Dlambda = log(10)*lambda; %  Derivative of lambda

n = size(Y, 1);

M = Z'*Z + lambda.*R;    %  perturbed cross-product matrix

Minv = inv(M);           %  inverse of this

ZMinv = Z*Minv;

A = ZMinv*Z';            %  smoothing or hat matrix

Yhat = A*Y;              %  matrix of smoothed values

Res  = Y - Yhat;         %  matrix of residuals

SSE = sum(sum(Res.^2));  %  error sum of squares

df = n - factor*trace(A);  %  degrees of freedom in the smooth

numGCV = SSE/n;          %  numerator of GCV

denGCV = (df/n)^2;       %  denominator of GCV

GCV = numGCV/denGCV;      %  GCV

DA = -Dlambda.*ZMinv*R*ZMinv'; %  derivative of A

DnumGCV = -2*trace(DA*Y*Res')/n;  %  derivative of SSE/n wrt lambda

DdenGCV = -2*df*factor*trace(DA)/n^2;

DGCV = ((denGCV*DnumGCV - numGCV*DdenGCV)/denGCV^2);

