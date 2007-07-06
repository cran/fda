function yinterp = safeinterp(x, y, xi, method, tol)
% SAVEINTERP uses INTERP1  for 1-D interpolation, but
%  protects it against values in xi being too close together.
% See help file for INTERP1 for further details.
% TOL is a nonnegative tolerance for closeness

%  Last modified 26 July 2006

if nargin < 5, tol = 1e-12;        end
if nargin < 4, method = 'linear';  end
if nargin < 3
    yinterp = y;
    return;
end

if tol < 0, error('TOL is negative.'); end

if isempty(diff(x) < tol)
    yinterp = interp1(x, y, xi, method);
else
    n = length(x);
    nonindex = [find(diff(x(:)) >= tol); n];
    xsafe = x(nonindex);
    ysafe = y(nonindex,:);
    yinterp = interp1(xsafe, ysafe, xi, method);
end

    