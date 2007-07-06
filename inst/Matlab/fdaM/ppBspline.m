function [Coeff,index] = ppBspline(t)
%PPBSPLINE computes the coefficients of the polynomials 
% of the piecewise polynomial B-spline of order k=length(t)-1
% corresponding to the knot sequence T 
% (i.e. B(k,T)(x) = (T(1+k)-T(1))[T(1),T(2),...,T(k+1)](.-x)_{+}^{k-1}),
% the coefficients of the polynomials each defines on 
% distinct intervals of T. 
% Say there are L distinct intervals in 
% T [T(knots[1]),T(knots[2])), [T(knots[2]),T(knots[3])),...,
% T(knots[L]),T(knots[L+1])], 
% then the coefficients are returned in the matrix COEFF as row vectors, 
% the i-th row corresponding to the coefficients of the polynomial on the
% interval [T(knots[i]),T(knots[i+1])), 
% such that 
% B(k,t)(x) = 
%   COEFF(i,1)*(x-T(i))^{k-1} + ... + COEFF(i,k-1)*(x-T(i)) + 
%   COEFF(i,k), for x in [T(knots[i]),T(knots[i+1])). 
% Note that we assume T(1) < T(k+1), i.e. T is not a
% sequence of the same knot. 
% The function returns the L times k matrix COEFF and the vector
% INDEX indicating where in the sequence T we find the knots. 
% This code uses part of the code from the function bspline.m .

%  Last modified 20 July 2006

k = length(t) - 1;
if k > 1
   adds  = ones(1,k-1);
   %  t(:) means t written as a column vector
   tt    = [adds*t(1), t(:)', adds*t(k+1)];
   a     = [adds*0,    1,     adds*0];
   inter = find(diff(tt)>0); 
   L     = length(inter);
   onesk = ones(1,2*(k-1));
   onesL = ones(L,1);
   tx    = onesL*(2-k:k-1) + inter'*onesk;
   tx(:) = tt(tx);
   tx    = tx - tt(inter)'*onesk;
   b     = onesL*(1-k:0)  + inter'*ones(1,k);
   b(:)  = a(b);
   %Coeff = sprpp(tx,b);
   km1 = k-1; 
   for r=1:km1
      for i=1:k-r
         b(:,i) =(tx(:,i+km1).*b(:,i)-tx(:,i+r-1).*b(:,i+1))./ ...
                 (tx(:,i+km1)-tx(:,i+r-1));
      end
   end
   Coeff = b;
   for r=2:k
      factor = (k-r+1)/(r-1);
      for i=k:-1:r
         Coeff(:,i) = (Coeff(:,i) - Coeff(:,i-1))*factor./tx(:,i+k-r);
      end
   end
   Coeff = Coeff(:,k:-1:1);
   index = inter - (k-1);
else
   Coeff = 1;
   index = 1;
end
