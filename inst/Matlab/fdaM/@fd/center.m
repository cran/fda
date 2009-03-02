function centerfd = center(fd)
%  CENTER Center functional observations by subtracting mean function
%  Returns CENTERFD, the centered functional data object

%  Last modified 2 March 2009

if ~isa_fd(fd)
    error ('Argument is not a functional data object.');
end

coef   = getcoef(fd);
coefd  = size(coef);
ndim   = length(coefd);
nrep   = coefd(2);
onebas = ones(1, nrep);
basisobj = getbasis(fd);
if ndim == 2
    coefmean = mean(coef,2);
    coef = coef - coefmean * onebas;
else
    nvar = coefd(3);
    for j = 1:nvar
        coefmean = mean(coef(:,:,j),2);
        coef(:,:,j) = coef(:,:,j) - coefmean * onebas;
    end
end
fdnames = getnames(fd);
if iscell(fdnames{3})
    fdnames{3}{1} = ['Centered ', fdnames{3}{1}];
else
    fdnames{3} = ['Centered ', fdnames{3}];
end

centerfd.coef     = coef;
centerfd.basisobj = basisobj;
centerfd.fdnames  = fdnames;

centerfd = class(centerfd, 'fd');

