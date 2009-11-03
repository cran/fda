function ijmat = ij2m(nrow,ncol)
%IJ2M sets up a NROW by NCOL matrix of indices of entries ordered row-wise
% example:  ij2m(2,3) produces
%            1  2  3
%            4  5  6
m = 0;
ijmat = zeros(nrow,ncol);
for i=1:nrow
    for j=1:ncol
        m = m + 1;
        ijmat(i,j) = m;
    end
end
