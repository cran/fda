function [wtvec, onewt, matwt] = wtcheck(n, wtvec)
%WTCHECK checks an input weight vector or weight matrix.  
%  If the weight vector is not empty, it checked for dimension, length and 
%  positivity.  In this case output ONEWT has value 0.
%  If the weight vector is empty, it is set up as a vector of ones of
%  the appropriate length.  In this case output ONEWT has value 1.
%  If weight is a matrix, it is checked for dimensions and being
%  positive definite, and matwt is set to 1, otherwise matwt is 0.

%  Last modified 9 July 2011 by Jim Ramsay

if nargin < 2,  wtvec = [];  end
if nargin < 1,  error('Argument n is missing.');  end

%  exchange arguments if wrong way around

if size(n(:),1) > 1 && size(wtvec(:),1) == 1
    tmp   = n;
    n     = wtvec;
    wtvec = tmp;
end

%  check n

if n ~= round(n)
    error('n is not an integer.');
end

if n < 1
    error('n is less than 1.');
end

%  check wtvec

if ~isempty(wtvec)
    sizew = size(wtvec);
    if any(isnan(wtvec(:)))
        error('WTVEC has NaN values.');
    end
    if all(sizew == n)
        %  WTVEC is a matrix of order n
        onewt = 0;
        matwt = 1;
        %  check weight matrix for being positive definite
        if issparse(wtvec)
            wteig = eig(full(wtvec));
        else
            wteig = eig(wtvec);
        end
        if any(strcmp(wteig, 'complex'))
            error('Weight matrix has complex eigenvalues.');
        end
        if min(wteig) <= 0
            error('Weight matrix is not positive definite.');
        end
    else
        %  WTVEC is treated as a vector
        if (length(sizew) > 1 && sizew(1) > 1 && sizew(2) > 1) || ...
                length(sizew) > 2
            error ('WTVEC is neither a vector nor a matrix of order n.');
        end
        wtvec = wtvec(:);
        if length(wtvec) == 1
            wtvec = wtvec.*ones(n,1);
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
