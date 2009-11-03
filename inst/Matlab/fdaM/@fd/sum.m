function sumfd = sum(fdobj)
%  SUM  Compute sum function for functional observations.
%  Argument:
%  FDOBJ ... a functional data object for the observations
%  Return:
%  SUMFD ... a functional data object for the sum

%  Last modified 20 July 2006

if ~isa_fd(fdobj)
    error ('Argument FDOBJ is not a functional data object.');
end

coef   = getcoef(fdobj);
coefd  = size(coef);
ndim   = length(coefd);
basisobj = getbasis(fdobj);
if ndim == 2
    coefsum = sum(coef,2);
else
    nvar = coefd(3);
    coefsum = zeros(coefd(1),1,nvar);
    for j = 1:nvar
        coefsum(:,1,j) = sum(coef(:,:,j),2);
    end
end

fdnames = getnames(fdobj);
fdnames{2} = 'Sum';
fdnames{3} = ['Sum ', fdnames{3}];
sumfd = fd(coefsum, basisobj, fdnames);

