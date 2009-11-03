function evalarray = eval_bifd(sevalarg, tevalarg, bifd, sLfd, tLfd)
%EVAL_BIFD  evaluates a bi-functional data object BIFD
%  at argument values in arrays SEVALARG and TEVALARG.
%  SLfd and TLfd are either integers giving the order of derivative,
%  or linear differential operators to be applied before evaluation.
%  Their defaults are 0, meaning that the function itself is evaluated.

%  last modified 20 January 2007

%  check that at least three arguments are present

if nargin < 3
    error('There are less than three arguments.');
end

%  exchange order if BIFD is the first argument

if isa_bifd(sevalarg)
    temp     = bifd;
    bifd     = sevalarg;
    sevalarg = tevalarg;
    tevalarg = temp;
end

%  check SEVALARG

sizesevalarg = size(sevalarg);
if sizesevalarg(1) > 1 && sizesevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
sevalarg = sevalarg(:);

%  check TEVALARG

sizetevalarg = size(tevalarg);
if sizetevalarg(1) > 1 && sizetevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
tevalarg = tevalarg(:);

if nargin < 5
    tLfd = int2Lfd(0);
end
if nargin < 4
    sLfd = int2Lfd(0);
end

sLfd = int2Lfd(sLfd);
tLfd = int2Lfd(tLfd);

if ~isa_bifd(bifd)
    error('Argument BIFD is not a bivariate functional data object.');
end
ns = length(sevalarg);
sbasisobj = getsbasis(bifd);
snbasis   = getnbasis(sbasisobj);

if ~isa_Lfd(sLfd)
    error ('sLfd is not a linear differential operator object.');
end

if ~isa_Lfd(tLfd)
    error ('sLfd is not a linear differential operator object.');
end

snderiv     = getnderiv(sLfd);

sbasismat = eval_basis(sevalarg, sbasisobj, snderiv);
if snderiv > 0 && ~strcmp(class(sLfd), 'double')
    derivwtmat = eval_basis(sevalarg, sbasisobj, sLfd);
    onerow = ones(1,snbasis);
    for j = 1:snderiv
        if any(abs(derivwtmat(:,j))) > 1e-7
            sbasismat = sbasismat +  ...
            (derivwtmat(:,j)*onerow) .* ...
                eval_basis(sevalarg, sbasisobj, j-1);
        end
    end
end

nt = length(tevalarg);
tbasisobj = gettbasis(bifd);
tnbasis   = getnbasis(tbasisobj);

tnderiv   = getnderiv(tLfd);

tbasismat = eval_basis(tevalarg, tbasisobj, tnderiv);
if tnderiv > 0 && ~strcmp(class(tLfd), 'double')
    derivwtmat = eval_basis(tevalarg, tbasisobj, tLfd);
    onerow = ones(1,tnbasis);
    for j = 1:tnderiv
        if any(abs(derivwtmat(:,j))) > 1e-7
            tbasismat = tbasismat +  ...
            (derivwtmat(:,j)*onerow) .* ...
                eval_basis(tevalarg, tbasisobj, j-1);
        end
    end
end

coef  = getcoef(bifd);
coefd = size(coef);
ndim  = length(coefd);

switch ndim
    case 2
        evalarray = sbasismat * coef * tbasismat';
    case 3
        ncurves  = coefd(3);
        evalarray = zeros(ns,nt,ncurves);
        for i= 1:ncurves
            evalarray(:,:,i) = sbasismat * coef(:,:,i) * tbasismat';
        end
    case 4
        ncurves  = coefd(3);
        nvar     = coefd(4);
        evalarray = zeros(ns,nt,ncurves,nvar);
        for i = 1:ncurves
            for j = 1:nvar
                evalarray(:,:,i,j) = sbasismat * coef(:,:,i,j) * tbasismat';
            end
        end
    otherwise
        error('coefficient array of improper dimension');
end
