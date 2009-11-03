function test = is_diag(c, EPS)
%  tests a matrix for being diagonal
%  C   ... matrix to be tested
%  EPS ... testing criterion: max off-diagonal element over min diagonal
%          element must be less than EPS

%  Last modified 2 December 2006

if nargin < 2
    EPS = 1e-12;
end
cd = size(c);
if length(cd) ~= 2
    test = 0;
else
    if cd(1) ~= cd(2)
        test = 0;
    else
        mindg  = min(abs(diag(c)));
        maxodg = max(max(abs(c - diag(diag(c)))));
        if mindg > 0
            if maxodg/mindg < EPS
                test = 1;
            else
                test = 0;
            end
        else
            if maxodg < EPS
                test = 1;
            else
                test = 0;
            end
        end
    end
end

