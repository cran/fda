function fdcell = fd2cell(fdobj)
%FD2CELL converts a univariate functional data object to a cell
%  object, mainly for purposes of defining a linear differential
%  operator object where the arguments are required to be cells.

%  Last modified 25 July 2006

%  check FDOBJ

if ~isa_fd(fdobj)
    error('FDOBJ is not a functional data object.');
end

%  get the coefficient matrix and the basis

coef     = getcoef(fdobj);
coefsize = size(coef);

%  check wether FDOBJ is univariate

if length(coefsize) > 2
    error('FDOBJ is not univariate.');
end
nrep     = coefsize(2);

fdcell = cell(1,nrep);
for i=1:nrep
    fdcell{i} = fdPar(fdobj(i));
end
