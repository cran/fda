function dy = derivs(tnow, y, wfd) 
% DERIVS sets up the 1st order system corresponding to   
%   linear differential operator defined by wfd.

% for example, if the system is of the second order,
%  y(1) is the current value of the variable
%  y(2) is the current value of the first derivative of the variable

%  last modified 31 October 2005

w  = eval_fd(tnow, wfd);
m  = length(w);
wmat = zeros(m, m);
wmat(1:(m-1),2:m) = eye(m-1);
wmat(m,:) = -w;
dy = wmat * y;

