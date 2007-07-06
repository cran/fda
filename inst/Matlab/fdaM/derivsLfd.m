function dy = derivsLfd(tnow, y) 
% DERIVS sets up the 1st order system corresponding to   
%   linear differential operator defined by wfd.

%  last modified 2 December 2006

global wfdcell;   

norder = length(wfdcell);
wmat   = zeros(norder, norder);
for j=1:(norder-1)
    wmat(j,j+1) = 1;
end

if norder == 1
    dy = eval_fd(tnow, -wfdcell{1})*y;
    return;
else
    for j=1:norder
        wj = eval_fd(tnow, wfdcell{j});
        wmat(norder,j) = -wj;
    end
    dy = wmat * y;
end
    

