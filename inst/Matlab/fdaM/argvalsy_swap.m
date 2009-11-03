function [y, argvals] = argvalsy_swap(argvals, y, basisobj)

%  carry out two tests, either of which calls for swapping
%  y and argvals

%  Range test:  Check ranges of both arguments, and
%  if y range is within range defined by basisobj, swap.

ymin = min(y(:));
ymax = max(y(:));

amin = min(argvals(:));
amax = max(argvals(:));

rangeval = getbasisrange(basisobj);

if (amin <  rangeval(1) || amax >  rangeval(2)) &&  ...
        (ymin >= rangeval(1) && ymax <= rangeval(2))
    ytmp    = y;
    y       = argvals;
    argvals = ytmp;
end

%  Dimension test:  If a3 has either 2nd or 3rd dimension length
%  greater than 1, and y3 has both 2nd and 3rd dimension lengths
%  equal to 1, swap.

if (size(y, 2)       == 1 && size(y, 2)       == 1) && ...
   (size(argvals, 2)  > 1 || size(argvals, 2)  > 1)
    ytmp    = y;
    y       = argvals;
    argvals = ytmp;
end


