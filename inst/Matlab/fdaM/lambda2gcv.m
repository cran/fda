function gcv = lambda2gcv(log10lambda, argvals, y, fdParobj, ...
                          wtvec, dffactor)
%  LAMBDA2GCV smooths data using smooth_basis and returns the 
%  GCV criterion that results.  DFFACTOR is an optional multiplier of df
%  in the definition of GCV

if nargin < 4, error('Insufficient number of arguments.');  end
if nargin < 6, dffactor = 1;               end
if nargin < 5, wtvec = ones(length(y),1);  end

fdParobj = putlambda(fdParobj, 10^log10lambda);

[fdobj, df, gcv] = ...
    smooth_basis(argvals, y, fdParobj, wtvec, dffactor);