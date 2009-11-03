function mmat = m2ij(nrow,ncol)
%M2IJ sets up a NROW*NCOL by 2 matrix of row-col indices associated
%  with a number of matrix entries row-wise
%  Example:  m2ij(2,3) produces
%     1     1
%     1     2
%     1     3
%     2     1
%     2     2
%     2     3
nval = nrow*ncol;
mmat = [reshape(ones(ncol,1)*(1:nrow), nval,1), ...
        reshape((1:ncol)'*ones(1,nrow),nval,1)];
