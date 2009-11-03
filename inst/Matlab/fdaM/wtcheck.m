function [wtvec, onewt, matwt] = wtcheck(n, wtvec)
%WTCHECK checks an input weight vector.  
%  If the weight vector is not empty, it checked for dimension, length and 
%  positivity.  In this case output ONEWT has value 0.
%  If the weight vector is empty, it is set up as a vector of ones of
%  the appropriate length.  In this case output ONEWT has value 1.

%  Last modified 8 November 2010 by Jim Ramsay

if nargin < 2,  wtvec = [];  end
if nargin < 1,  error('Argument n is missing.');  end

if n ~= round(n)
    error('n is not an integer.');
end

if n < 1
    error('n is less than 1.');
end

if ~isempty(wtvec)
    sizew = size(wtvec);
    if all(sizew == n)
        %  WTVEC is a matrix of order n
        onewt = 0;
        matwt = 1;
    else
        %  WTVEC is treated as a vector
        if (length(sizew) > 1 && sizew(1) > 1 && sizew(2) > 1) || ...
                length(sizew) > 2
            error ('WTVEC must be a vector.');
        end
        if length(sizew) == 2 && sizew(1) == 1
            wtvec = wtvec';
        end
        if length(wtvec) ~= n
            error('WTVEC of wrong length');
        end
        if min(wtvec) <= 0
            error('All values of WTVEC must be positive.');
        end
        onewt = 0;
        matwt = 0;
    end
else
    wtvec = ones(n,1);
    onewt = 1;
    matwt = 0;
end
