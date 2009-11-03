function [argvals, n] = argcheck(argvals)

%  check ARGVALS

if ~strcmp(class(argvals), 'double')
    error('ARGVALS is not of class double.');
end

if size(argvals,1) == 1
    argvals = argvals';
end

[n, ncl] = size(argvals);  %  number of observations

if ncl > 1
    error('ARGVALS is not a vector.')
end
if n < 2
    error('ARGVALS does not contain at least two values.');
end

