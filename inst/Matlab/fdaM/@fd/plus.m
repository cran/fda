function plusfd = plus(fd1, fd2, basisobj)
%  PLUS: Pointwise sum of two functional data objects,
%    the sum of a scalar and a functional data object,
%    or the sum of a vector and a functional data obect
%       where the length of the vector is the same as the
%       number of replications of the object.
%  When both arguments are functional data objects, 
%  they need not have the same bases, 
%  but they must either (1)  have the same number of replicates, or
%  (2) one function must have a single replicate and other multiple 
%  replicates.  In the second case, the singleton function is 
%  replicated to match the number of replicates of the other function.
%  In either case, they must have the same number of functions.
%  When both arguments are functional data objects, and the
%  bases are not the same,
%  the basis used for the sum is constructed to be of higher 
%  dimension than the basis for either factor according to rules
%  described in function TIMES for two basis objects. 
%  Finally, in the simple case where both arguments are
%  functional data objects, the bases are the same, and the
%  coefficient matrices are the same sizes, the coefficient
%  matrices are simply added.

%  last modified 1 December 2006

if ~(isa_fd(fd1) || isa_fd(fd2))
      error('Neither argument for + is a functional data object.');
end

if isa_fd(fd1) && isa_fd(fd2)
    %  both arguments are functional data objects
    %  check to see of the two bases are identical
    %  and if the coefficient matrices are conformable.
    basisobj1 = getbasis(fd1);
    basisobj2 = getbasis(fd2);
    type1   = getbasistype(basisobj1);
    type2   = getbasistype(basisobj2);
    nbasis1 = getnbasis(basisobj1);
    nbasis2 = getnbasis(basisobj2);
    range1  = getbasisrange(basisobj1);
    range2  = getbasisrange(basisobj2);
    params1 = getbasispar(basisobj1);
    params2 = getbasispar(basisobj2);
    coef1  = getcoef(fd1);
    coef2  = getcoef(fd2);
    coefd1 = size(coef1);
    coefd2 = size(coef2);
    %  test to see if the two objects match completely
    if strcmp(type1, type2)         && ...
            all(range1  == range2)  && ...
            nbasis1 == nbasis2      && ...
            all(params1 == params2) && ...
            all(coefd1 == coefd2)
        %  the two coefficient matrices can be simply added
        fdnames = getnames(fd1);
        plusfd.coef = coef1 + coef2;
        plusfd.basisobj = basisobj1;
        plusfd.fdnames  = fdnames;
        plusfd = class(plusfd, 'fd');
        return;
    end
    %  check to see if the number of dimensions match
    ndim1  = length(coefd1);
    ndim2  = length(coefd2);
    if ndim1 ~= ndim2
        error('Dimensions of coefficient matrices not compatible.');
    end
    %  allow for one function being a single replicate,
    %  and if so, copy it as many times as there are replicates
    %  in the other function.
    if coefd1(2) == 1 && coefd2(2) >  1
        if     ndim1 == 2
            coef1 = coef1*ones(1,coefd2(2));
        elseif ndim1 == 3
            temp = zeros(coefd2);
            for j=1:coefd1(3)
                temp(:,:,j) = squeeze(coef1(:,1,j))*ones(1,coefd2(2));
            end
            coef1 = temp;
        else
            error('Dimensions of coefficient matrices not compatible.');
        end
        coefd1 = size(coef1);
        fd1    = putcoef(fd1, coef1);
    end
    if coefd1(2) >  1 && coefd2(2) == 1
        if     ndim2 == 2
            coef2 = coef2*ones(1,coefd1(2));
        elseif ndim1 == 3
            temp = zeros(coefd1);
            for j=1:coefd2(3)
                temp(:,:,j) = squeeze(coef2(:,1,j))*ones(1,coefd1(2));
            end
            coef2 = temp;
        else
            error('Dimensions of coefficient matrices not compatible.');
        end
        coefd2 = size(coef2);
        fd2    = putcoef(fd2, coef2);
    end
    %  check for equality of dimensions of coefficient matrices
    if coefd1(2) ~= coefd2(2) 
        error('Number of replications are not equal.');
    end
    %  check for equality of numbers of functions
    if ndim1 > 2 && ndim2 > 2 && ndim1 ~= ndim2
        error(['Both arguments multivariate, ',  ...
               'but involve different numbers ', ...
               'of functions.']);
    end
    basisobj1 = getbasis(fd1);
    basisobj2 = getbasis(fd2);
    %  check for equality of two bases
    if basisobj1 == basisobj2
        %  if equal, just difference coefficient matrices
        fdnames = getnames(fd1);
        plusfd.coef = coef1 + coef2;
        plusfd.basisobj = basisobj1;
        plusfd.fdnames  = fdnames;
        plusfd = class(plusfd, 'fd');
        return;
    else
        nbasis1   = getnbasis(basisobj1);
        nbasis2   = getnbasis(basisobj2);
        rangeval1 = getbasisrange(basisobj1);
        rangeval2 = getbasisrange(basisobj2);
        if (any(rangeval1 ~= rangeval2))
            error('The ranges of the arguments are not equal.');
        end
        neval     = max(10*max(nbasis1+nbasis2) + 1, 101);
        evalarg   = linspace(rangeval1(1), rangeval2(2), neval)';
        fdarray1  = eval_fd(fd1, evalarg);
        fdarray2  = eval_fd(fd2, evalarg);
        if (ndim1 <= 2 && ndim2 <= 2) || (ndim1 > 2 && ndim2 > 2)
            fdarray = fdarray1 + fdarray2;
        end
        if ndim1 == 2 && ndim2 > 2
            fdarray = zeros(coefd2);
            for ivar = 1:coefd2(3)
                fdarray(:,:,ivar) = fdarray1 + fdarray2(:,:,ivar);
            end
        end
        if ndim1 > 2 && ndim2 == 2
            fdarray = zeros(coefd1);
            for ivar = 1:coefd1(3)
                fdarray(:,:,ivar) = fdarray1(:,:,ivar) + fdarray2;
            end
        end
        %  set up basis for sum
        if nargin < 3, basisobj = basisobj1.*basisobj2;  end
        coefsum = project_basis(fdarray, evalarg, basisobj, 1);
        fdnames1 = fd1.fdnames;
        fdnames2 = fd2.fdnames;
        fdnames  = fdnames1;
        fdnames{3} = [fdnames1{3},'+',fdnames2{3}];
    end
 else
    %  one argument is numeric and the other is functional
    if ~(isnumeric(fd1) || isnumeric(fd2))
        error('Neither argument for + is numeric.');
    end
    if isnumeric(fd1) && isa_fd(fd2)
        fac = fd1;
        fd  = fd2;
    elseif isa_fd(fd1) && isnumeric(fd2)
        fac = fd2;
        fd  = fd1;
    else
        error('One of the arguments for + is of the wrong class.');
    end
    coef     = getcoef(fd);
    coefd    = size(coef);
    basisobj = getbasis(fd);
    nbasis   = getnbasis(basisobj);
    rangeval = getbasisrange(basisobj);
    neval    = max(10*nbasis + 1,201);
    neval    = min(neval,201);
    evalarg  = linspace(rangeval(1),rangeval(2), neval)';
    fdmat    = eval_fd(evalarg, fd);
    %  If one of the objects has length 1 and the other
    %  is longer, expand the scalar object into a vector
    if length(fac) ~= coefd(2)
        if length(fac) == 1 && coefd(2) > 1
            % nothing needs to be done in this case
        elseif length(fac) > 1 && coefd(2) == 1
            fdmat = fdmat*ones(1,length(fac));
            fac   = ones(neval,1)*fac(:)';
        else
            error(['Dimensions of numerical factor and functional', ...
                   ' factor cannot be reconciled.']);
        end
    end
    fdarray  = fac + fdmat;
    coefsum = project_basis(fdarray, evalarg, basisobj);
    fdnames  = fd.fdnames;
    if length(fac) == 1
        fdnames{3} = [num2str(fac),' + ',fdnames{3}];
    end
end

plusfd.coef     = coefsum;
plusfd.basisobj = basisobj;
plusfd.fdnames  = fdnames;

plusfd = class(plusfd, 'fd');