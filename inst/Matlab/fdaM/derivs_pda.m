function dy = derivs(tnow, y, bwtcell, awtcell, ufdcell) 
% DERIVS sets up the 1st order system corresponding to   
%   linear differential operator defined by wfd.

% for example, if the system is of the second order,
%  y(1) is the current value of the variable
%  y(2) is the current value of the first derivative of the variable

%  last modified 25 January 2011

if nargin < 5, ufdcell = [];  end
if nargin < 4, awtcell = [];  end

m    = length(bwtcell);
wmat = zeros(m, m);
wmat(1:(m-1),2:m) = eye(m-1);
for j=1:m
    wmat(m,j) = -eval_fd(tnow, getfd(bwtcell{j}));
end
dy = wmat * y;
if ~isempty(awtcell) && ~isempty(ufdcell)
    n = length(ufdcell);
    for k=1:n
        dy(m) = dy(m) + eval_fd(tnow, getfd(awtcell{k}));
    end
end

