function  [res, U, indlo, indhi] = ...
                              reschk(res, rangey, sigma0) 
% RESCHK brings residuals outside of limits to limits
%  First center the residuals.  See discussion of zmat
%  about the conditions that zmat must satisfy to make this
%  operation legitimate.
res = res - mean(res);
%  Look for residuals below lower boundary
indlo = find(res < rangey(1)*sigma0);
if length(indlo) > 0
    res(indlo) = rangey(1)*sigma0;
end
%  Look for residuals above upper boundary
indhi = find(res > rangey(2)*sigma0);
if length(indhi) > 0
    res(indhi) = rangey(2)*sigma0;
end
%  Compute standardized residuals
U   = res./sigma0;
