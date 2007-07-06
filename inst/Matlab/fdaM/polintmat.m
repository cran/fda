function [y, dy] = polintmat(xa, ya, x)
%  YA is an 3-D array with 1st D same as XA
n  = length(xa);
ya = squeeze(ya);
yad = size(ya);
if yad(1) ~= n
    error('First dimension of YA must match XA');
end
difx = xa - x;
absxmxa = abs(difx);
nseq = 1:n;
ns = min(nseq(absxmxa == min(absxmxa)));
ds = ya;
cs = ya;
y  = ya(ns,:);
ns = ns - 1;
for m = 1:(n-1)
    for i = 1:(n-m)
        ho      = difx(i);
        hp      = difx(i+m);
        w       = (cs(i+1,:) - ds(i,:))./(ho - hp);
        ds(i,:) = hp.*w;
        cs(i,:) = ho.*w;
    end
    if 2*ns < n-m
        dy = cs(ns+1,:);
    else
        dy = ds(ns,:);
        ns = ns - 1;
    end
    y = y + dy;
end

