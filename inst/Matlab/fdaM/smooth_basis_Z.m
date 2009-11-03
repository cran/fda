function [fdobj, beta, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_basis_Z(argvals, y, fdParobj, Z, wtvec, dffactor, fdnames)
%SMOOTH_BASIS_Z smooths data by penalized least squares.  The smooth curves
%  are expressed as a basis function expansions, and this function 
%  computes the coefficients of the expansions.  Smoothness is controlled
%  by controlling the size of an integrated squared linear differential
%  operator.  The integral is multiplied by a smoothing or bandwidth
%  parameter.
%  In addition, a matrix of covariate values may be input to augment the
%  smooth component of the fit.  A smooth augmented in this way is often
%  called a semiparametric smooth.
%
%  Required arguments for this function are:
%
%  ARGVALS  ... A set of argument values, set by default to equally spaced
%               on the unit interval (0,1).
%  Y        ... an array containing values of curves
%               If the array is a matrix, rows must correspond to argument
%               values and columns to replications, and it will be assumed
%               that there is only one variable per observation.
%               If Y is a three-dimensional array, the first dimension
%               corresponds to argument values, the second to replications,
%               and the third to variables within replications.
%               If Y is a vector, only one replicate and variable are 
%               assumed.
%  FDPAROBJ ... A functional parameter or fdPar object.  This object 
%               contains the specifications for the functional data
%               object to be estimated by smoothing the data.  See
%               comment lines in function fdPar for details.
%               This argument may also be either a FD object, or a 
%               BASIS object.  If this argument is a basis object, the 
%               smoothing parameter LAMBDA is set to 0.
%
%  Optional arguments are:
%
%  Z            A design matrix for semiparametric regression.
%               It must be n by q, where n is length(argvals).
%  WTVEC    ... A vector of N weights, set to one by default, that can
%               be used to differentially weight observations.
%  DFFACTOR ... A multiplier of df in GCV, set to one by default
%  FDNAMES  ... A cell of length 3 with names for
%               1. argument domain, such as 'Time'
%               2. replications or cases
%               3. the function.
%  Returned objects are:
%    FDOBJ  ...  an object of class fd containing coefficients.
%    BETA   ...  the regression coefficients for the covariates if supplied
%                or empty otherwise
%    DF     ...  a degrees of freedom measure.
%    GCV    ...  a measure of lack of fit discounted for df.
%                If the function is univariate, GCV is a vector 
%                containing the error  sum of squares for each 
%                function, and if the function is multivariate, 
%                GCV is a NVAR by NCURVES matrix.
%    COEF   ...  The coefficient matrix for the basis function
%                  expansion of the smoothing function in the 
%                first NBASIS columns.
%                If a design matrix is included with q columns,
%                the corresponding regression coefficients are
%                in the last q columns of COEF.
%    SSE    ...  the error sums of squares.  
%                SSE is a vector or matrix of the same size as 
%                GCV.
%    PENMAT ...  the penalty matrix.
%    Y2CMAP ...  the matrix mapping the data to the coefficients.

%  Last modified 20 July 2010 by Jim Ramsay

if nargin < 3
    error('There is not at least three arguments.');
end

%  check ARGVALS

if ~strcmp(class(argvals), 'double')
    error('ARGVALS is not of class double.');
end

if size(argvals,1) == 1
    argvals = argvals';
end

[n, ncl] = size(argvals);  %  number of observations
if ncl > 1
    error('ARGVALS is not a vector.')
end
if n < 2
    error('ARGVALS does not contain at least two values.');
end

%  check Y

if ~strcmp(class(y), 'double')
    error('Y is not of class double.');
end

ydim = size(y);
if length(ydim) == 2 && ydim(1) == 1
    y = y';
end

ydim = size(y);  %  number of observations
if ydim(1) ~= n
    error('Y is not the same length as ARGVALS.');
end

%  set default argument values

if nargin < 7
    fdnames{1} = 'arguments';
    fdnames{2} = 'replications';
    fdnames{3} = 'variables';
end

if nargin < 6, dffactor = 1;         end

if nargin < 5, wtvec  = ones(n,1);   end

if nargin < 4, Z = [];               end

%  check fdParobj

if ~isa_fdPar(fdParobj) 
    if isa_fd(fdParobj) || isa_basis(fdParobj)
        fdParobj = fdPar(fdParobj);
    else
        error(['FDPAROBJ is not a functional parameter object, ', ...
               'not a functional data object, and ', ...
               'not a basis object.']);
    end
end

%  check LFDOBJ

Lfdobj = getLfd(fdParobj);
Lfdobj = int2Lfd(Lfdobj);

%  check BASIS

fdobj    = getfd(fdParobj);
basisobj = getbasis(fdobj);
if ~isa_basis(basisobj)
    error('BASIS is not a basis object.');
end

nbasis   = getnbasis(basisobj) - length(getdropind(basisobj));

%  check design matrix Z

if ~isempty(Z)
    [nZ, q] = size(Z);
    if nZ ~= n
        error('Design matrix Z does not have n rows.');
    end
else
    q = 0;
end

%  check WTVEC

if ~isempty(wtvec)
    sizew = size(wtvec);
    if (length(sizew) > 1 && sizew(1) > 1 && sizew(2) > 1) || ...
            length(sizew) > 2
        error ('WTVEC must be a vector.');
    end
    if length(sizew) == 2 && sizew(1) == 1
        wtvec = wtvec';
    end
    if length(wtvec) ~= n
        error('WTVEC of wrong length');
    end
    if min(wtvec) <= 0
        error('All values of WTVEC must be positive.');
    end
else
    wtvec = ones(n,1);
end

if isempty(dffactor), dffactor = 1.0;  end

%  check LAMBDA

lambda = getlambda(fdParobj);
if lambda < 0
    warning ('Wid1:negative',...
        'Value of LAMBDA was negative, and 0 used instead.');
    lambda = 0;
end

%  -----------------------------------------------------------------------
%                      Set up analysis
%  -----------------------------------------------------------------------

%  set number of curves and number of variables

sizey = size(y);
ndim  = length(sizey);
switch ndim
    case 1
        ncurves = 1;
        nvar    = 1;
    case 2
        ncurves = sizey(2);
        nvar    = 1;
    case 3
        ncurves = sizey(2);
        nvar    = sizey(3);
    otherwise
        error('Second argument must not have more than 3 dimensions');
end

%  set up matrix of basis function values

basismat  = eval_basis(argvals, basisobj);

%  ------------------------------------------------------------------
%                set up the linear equations for smoothing
%  ------------------------------------------------------------------

if n >= nbasis || lambda > 0
    
    %  The following code is for the coefficients completely determined
    
    BZmat  = basismat;
    if q > 0
        BZmat = [BZmat, Z];
    end
    BZmatw    = BZmat .* (wtvec * ones(1,nbasis+q));
    BZtBZmat  = BZmatw' * BZmat;
    BZtBZmat0 = BZtBZmat;
    
    %  set up right side of equation
    
    if ndim < 3
        Dmat = BZmatw' * y;
    else
        Dmat = zeros(nbasis,ncurves,nvar);
        for ivar = 1:nvar
            Dmat(:,:,ivar) = BZmatw' * y(:,:,ivar);
        end
    end
    
    %  set up regularized cross-product matrix BMAT

    if lambda > 0
        penmat = eval_penalty(basisobj, Lfdobj);
        Bnorm   = sqrt(sum(sum(Bmat.^2)));
        pennorm = sqrt(sum(sum(penmat.^2)));
        condno  = pennorm/Bnorm;
        if lambda*condno > 1e12
            lambda = 1e12/condno;
            warning('Wid2:reduce', ...
                ['LAMBDA reduced to ',num2str(lambda), ...
                    ' to prevent overflow']);
        end
        if q == 0
            BZtBZmat = BZtBZmat + lambda .* penmat;
        else
            ind = 1:nbasis;
            BZtBZmat(ind,ind) = BZtBZmat(ind,ind) + lambda .* penmat;
        end
    end
    
    %  compute inverse of Bmat
    
    if is_diag(BZtBZmat)
        BZtBZmatinv = diag(1./diag(BZtBZmat));
    else
        BZtBZmatinv = inv(BZtBZmat);
    end
    
    %  ------------------------------------------------------------------
    %       Compute the coefficients defining the smooth and 
    %            summary properties of the smooth
    %  ------------------------------------------------------------------

    %  compute map from y to c
    
    y2cMap = BZtBZmatinv * BZtBZw';
    
    %  compute degrees of freedom of smooth
    
    df = full(sum(diag(BZtBZmatinv * BZtBZmat0)));
    
    %  solve normal equations for each observation
    
    if ndim < 3
        coef = BZtBZmatinv * BZtBZmat;
    else
        coef = zeros(nbasis+q, ncurves, nvar);
        for ivar = 1:nvar
            coef(:,:,ivar) = BZtBZmatinv * Dmat(:,:,ivar);
        end
    end
    
else
    
    error('Solution is under-determined.');
    
end

%  ------------------------------------------------------------------
%            compute SSE, yhat, GCV and other fit summaries
%  ------------------------------------------------------------------

%  compute error sum of squares

if ndim < 3
    yhat = BZtBZmat * coef;
    SSE = sum((y - yhat).^2);
else
    SSE = zeros(nvar,ncurves);
    for ivar = 1:nvar
        coefi = squeeze(coef(:,:,ivar));
        yhati = BZtBZmat * coefi;
        yi    = squeeze(y(:,:,ivar));
        SSE(ivar,:) = sum((yi - yhati).^2);
    end
end

%  compute  GCV index

if df < n
    gcv = (SSE./n)./((n - dffactor*df)/n)^2;
else
    gcv = NaN;
end

if q == 0
    fdobj = fd(coef, basisobj, fdnames);
else
    ind   = 1:nbasis;
    if ndim < 3
        fdobj = fd(coef(ind,:), basisobj, fdnames);
    else
        fdobj = fd(coef(ind,:,:), basisobj, fdnames);
    end
end

%  set up the regression coefficient matrix beta

if q > 0
    ind = (nbasis+1):(nbasis+q);
    if ndim < 3
        beta = coef(ind,:);
    else
        beta = coef(ind,:,:);
    end
else
    beta = [];
end



