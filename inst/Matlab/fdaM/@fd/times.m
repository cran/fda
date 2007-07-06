function timesfd = times(fdobj1, fdobj2, basisobj)
%  TIMES: Pointwise product of two functional data objects,
%    the product of a scalar and a functional data object,
%    or the product of a vector and a functional data obect
%       where the length of the vector is the same as the
%       number of replications of the object.
%  When both arguments are functional data objects, 
%  they need not have the same bases, 
%  but they must either (1)  have the same number of replicates, or
%  (2) one function must have a single replicate and other multiple 
%  replicates.  In the second case, each function in the multiple
%  replicate object is multiplied by the singleton function in the 
%  other objects.
%  In either case, they must have the same number of functions.

%  When both arguments are functional data objects, the
%  basis used for the product is constructed to be of higher 
%  dimension than the basis for either factor according to rules
%  described in function TIMES for two basis objects.  

%  last modified 12 December 2006

%  Check if at least one argument is a functional data object

if ~(isa_fd(fdobj1) || isa_fd(fdobj2))
      error('Neither argument for * is a functional data object.');
end

%  Determine which of two cases hold:
%   1.  both variables are functional
%   2.  only one of them is functional

if isa_fd(fdobj1) && isa_fd(fdobj2)
    
    %  --------------------------------------------------------
    %       both arguments are functional data objects
    %  --------------------------------------------------------
    
    %  get the dimensions of the two objects
    
    coef1  = getcoef(fdobj1);
    coef2  = getcoef(fdobj2);
    coefd1 = size(coef1);
    coefd2 = size(coef2);
    ndim1  = length(coefd1);
    ndim2  = length(coefd2);
    
    %  check that the two coefficient arrays have the same
    %  number of dimensions
    
    if ndim1 ~= ndim2
        error('Dimensions of coefficient matrices not compatible.');
    end
    
    %  allow for one function having a single replicate,
    %  and if so, copy it as many times as there are replicates
    %  in the other function.
    
    %  FDOBJ1 is single,  FDOBJ2 has replications
    
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
        fdobj1    = putcoef(fdobj1, coef1);
    end
    
    %  FDOBJ2 is single,  FDOBJ1 has replications
    
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
        fdobj2    = putcoef(fdobj2, coef2);
    end
    
    %  Now check that the numbers of replications match
    
    if coefd1(2) ~= coefd2(2) 
        error('Number of replications are not equal.');
    end
    
    %  check for matching in the multivariate case
    
    if ndim1 > 2 && ndim2 > 2 && ndim1 ~= ndim2
        error(['Both arguments multivariate, ',  ...
                'but involve different numbers ', ...
                'of functions.']);
    end
    
    %  extract the two bases
    
    basisobj1 = getbasis(fdobj1);
    basisobj2 = getbasis(fdobj2);
    
    %  check that the ranges match if a range not supplied
    
    rangeval1 = getbasisrange(basisobj1);
    rangeval2 = getbasisrange(basisobj2);
    if (any(rangeval1 ~= rangeval2))
        error('The ranges of the arguments are not equal.');
    end
    rng = rangeval1;
    
    %  set default basis object
    
    if nargin < 3
        basisobj = basisobj1.*basisobj2; 
    end
    
    %  set up a fine mesh for evaluating the product
    
    neval   = max(10*getnbasis(basisobj), 101);
    evalarg = linspace(rng(1), rng(2), neval)';
    
    %  set up arrays of function values
    
    fdarray1 = eval_fd(evalarg, fdobj1);
    fdarray2 = eval_fd(evalarg, fdobj2);
    
    %  compute product arrays
    
    if (ndim1 <= 2 && ndim2 <= 2) || (ndim1 > 2 && ndim2 > 2)
        %  product array where the number of dimensions match
        fdarray = fdarray1.*fdarray2;
    else
        %  product array where the number of dimensions don't match
        if ndim1 == 2 && ndim2 > 2
            fdarray = zeros(coefd2);
            for ivar = 1:coefd2(3)
                fdarray(:,:,ivar) = ...
                    fdarray1 .* squeeze(fdarray2(:,:,ivar));
            end
        end
        if ndim1 > 2 && ndim2 == 2
            fdarray = zeros(coefd1);
            for ivar = 1:coefd1(3)
                fdarray(:,:,ivar) = ...
                    squeeze(fdarray1(:,:,ivar)) .* fdarray2;
            end
        end
    end
    
    %  set up the coefficient by projecting on to the 
    %  product basis
    
    coefprod = project_basis(fdarray, evalarg, basisobj, 1);
    
    %  set up the names
    
    fdnames1 = fdobj1.fdnames;
    fdnames2 = fdobj2.fdnames;
    fdnames  = fdnames1;
    fdnames{3} = [fdnames1{3},'*',fdnames2{3}];
    
else
    
    %  --------------------------------------------------------
    %    one argument is numeric and the other is functional
    %  --------------------------------------------------------
    
    if ~(isnumeric(fdobj1) || isnumeric(fdobj2))
        error('Neither argument for * is numeric.');
    end
    if isnumeric(fdobj1) && isa_fd(fdobj2)
        fac   = fdobj1;
        fdobj = fdobj2;
    elseif isa_fd(fdobj1) && isnumeric(fdobj2)
        fac   = fdobj2;
        fdobj = fdobj1;
    else
        error('One of the arguments for .* is of the wrong class.');
    end
    coef  = getcoef(fdobj);
    if length(size(coef)) == 2
        [nbasis, N] = size(coef);
        nvar = 1;
    else
        [nbasis, N, nvar] = size(coef);
    end
    if ~(length(fac) == N || length(fac) == 1)
        error('The length of the numerical factor is incorrect.');
    end
    fdnames  = fdobj.fdnames;
    if length(fac) == 1
        fdnames{3} = [num2str(fac),'*',fdnames{3}];
        fac = fac.*ones(N,1);
    end
    [n1,n2] = size(fac);
    if ~(n1 == 1 || n2 == 1)
        error('The numerical factor is not a vector.');
    end
    if n1 ==1
        fac = fac';
    end
    facmat = ones(nbasis,1)*fac';
    if length(size(coef)) == 2
        coefprod = facmat.*coef;
    else
        coefprod = zeros(nbasis, N, nvar);
        for j=1:nvar
            coefprod(:,:,j) = facmat.*squeeze(coef(:,:,j));
        end
    end
    basisobj = getbasis(fdobj);
end

timesfd.coef     = coefprod;
timesfd.basisobj = basisobj;
timesfd.fdnames  = fdnames;

timesfd = class(timesfd, 'fd');