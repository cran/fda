function basisobj = create_basis(type, rangeval, nbasis, params, ...
                                 dropind)
%  CREATE_BASIS  An alternative call to CREATE_BASIS_fd.

if nargin < 5, dropind = [];  end

basisobj = basis(type, rangeval, nbasis, params, dropind);


