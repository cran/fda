function evalarray = evaldiag_bifd(evalarg, bifdobj, sLfd, tLfd)
%EVALDIAG_BIFD  evaluates a bi-functional data object BIFD
%  with both argument values in array EVALARG.
%  SLfd and TLfd are either integers giving the order of derivative,
%  or linear differential operators to be applied before evaluation.
%  Their defaults are 0, meaning that the function itself is evaluated.

%  last modified 20 July 2006

%  check that at least three arguments are present

if nargin < 2
    error('There are less than two arguments.');
end

%  exchange order if BIFD is the first argument

if isa_bifd(evalarg)
    temp    = bifdobj;
    bifdobj = evalarg;
    evalarg = temp;
end

%  check EVALARG

sizeevalarg = size(evalarg);
if sizeevalarg(1) > 1 && sizeevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
evalarg = evalarg(:);

%  set default differential operators

if nargin < 4
    tLfd = int2Lfd(0);
end

if nargin < 3
    sLfd = int2Lfd(0);
end

if ~isa_bifd(bifdobj)
    error('Argument BIFD is not a bivariate functional data object.');
end

n = length(evalarg);

%  extract the two bases

sbasisobj = getsbasis(bifdobj);
tbasisobj = gettbasis(bifdobj);

%  check that the bases have the same range

ranges = getbasisrange(sbasisobj);
ranget = getbasisrange(tbasisobj);
if any(ranges ~= ranget)
    error('The ranges are not identical.');
end

%  check the differential operators

sLfd = int2Lfd(sLfd);
tLfd = int2Lfd(tLfd);

if ~isa_Lfd(sLfd)
    error ('SLFD is not a linear differential operator object.');
end

if ~isa_Lfd(tLfd)
    error ('TLFD is not a linear differential operator object.');
end

%  compute the basis matrix for SBASISOBJ

snderiv   = getnderiv(sLfd);
if snderiv > 0 && ~strcmp(class(sLfd), 'double')
    derivwtmat = eval_basis(evalarg, sbasisobj, sLfd);
    for j = 1:snderiv
        if any(abs(derivwtmat(:,j))) > 1e-7
            sbasismat <- sbasismat +  ...
            (derivwtmat(:,j)*onerow) .* ...
                eval_basis(evalarg, sbasisobj, j-1);
        end
    end
end

%  compute the basis matrix for tBASISOBJ

tnderiv   = getnderiv(tLfd);
if tnderiv > 0 && ~strcmp(class(tLfd), 'double')
    derivwtmat = eval_basis(evalarg, tbasisobj, tLfd);
    for j = 1:tnderiv
        if any(abs(derivwtmat(:,j))) > 1e-7
            tbasismat <- tbasismat +  ...
            (derivwtmat(:,j)*onerow) .* ...
                eval_basis(evalarg, tbasisobj, j-1);
        end
    end
end

%  Extract the coefficient matrix

coef  = getcoef(bifdobj);
coefd = size(coef);
ndim  = length(coefd);

switch ndim
    case 2
        evalarray = diag(sbasismat * coef * tbasismat');
    case 3
        ncurves   = coefd(3);
        evalarray = zeros(n,ncurves);
        for i=1:ncurves
            evalarray(:,i) = diag(sbasismat * coef(:,:,i) * tbasismat');
        end
    case 4
        ncurves  = coefd(3);
        nvar     = coefd(4);
        evalarray = zeros(n,ncurves,nvar);
        for i = 1:ncurves
            for j = 1:nvar
                evalarray(:,i,j) = ...
                    diag(sbasismat * coef(:,:,i,j) * tbasismat');
            end
        end
    otherwise
        error('coefficient array of improper dimension');
end
