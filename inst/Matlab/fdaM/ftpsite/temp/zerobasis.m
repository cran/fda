function zerobasmat = zerobasis(k)
%  sets up a k by k-1 matrix with orthonormal columns
%  using the first k non-constant fourier basis function
%  values at 0.5, ..., k-0.5
fbasis = create_fourier_basis(k,k);
tk = (0:k-1) + 0.5;
fbasmat = eval_basis(tk, fbasis);
fbasmat = fbasmat(:,2:k);
fbasnorm = sum(fbasmat.^2);
zerobasmat = fbasmat./(ones(k,1)*fbasnorm);



