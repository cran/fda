function geigenstr = geigen(Amat, Bmat, Cmat)
%  GEIGEN  solve the generalized eigenanalysis problem
%
%    max tr L'AM / sqrt(tr L'BL tr M'CM) w.r.t. L and M
%
%  Arguments:
%  AMAT ... p by q matrix
%  BMAT ... order p symmetric positive definite matrix
%  CMAT ... order q symmetric positive definite matrix
%  Returns:
%  VALUES ... vector of length s = min(p,q) of eigenvalues
%  LMAT   ... p by s matrix L
%  MMAT   ... q by s matrix M

%  last modified 20 July 2006

Bsize = size(Bmat);
Csize = size(Cmat);
if Bsize(1) ~= Bsize(2)
    error('BMAT is not square');
end
if (Csize(1) ~= Csize(2))
    error('CMAT is not square');
end
p = Bsize(1);
q = Csize(1);
if max(max(abs(Bmat - Bmat')))/max(max(abs(Bmat))) > 1e-10
    error('BMAT not symmetric.');
else
    Bmat = (Bmat + Bmat')./2;
end
if max(max(abs(Cmat - Cmat')))/max(max(abs(Cmat))) > 1e-10
    error('CMAT not symmetric.');
else
    Cmat = (Cmat + Cmat')./2;
end
Bfac  = chol(Bmat);
Cfac  = chol(Cmat);
Bfacinv = inv(Bfac);
Cfacinv = inv(Cfac);
Dmat = Bfacinv' * Amat * Cfacinv;
if (p >= q)
    [u,d,v] = svd(Dmat);
    values = d;
    Lmat = Bfacinv * u;
    Mmat = Cfacinv * v;
else
    [u,d,v] = svd(Dmat');
    values = d;
    Lmat = Bfacinv * v;
    Mmat = Cfacinv * u;
end

geigenstr.values = values;
geigenstr.Lmat   = Lmat;
geigenstr.Mmat   = Mmat;

