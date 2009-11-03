function [fdobj, df, gcv, SSE, penmat, y2cMap, argvals, y] = ...
    smooth_basis(argvals, y, fdParobj, wtvec, fdnames, dffactor)
%SMOOTH_BASIS  Smooths data by penalized least squares.  The smooth curves
%  are expressed as a basis function expansions, and this function 
%  computes the coefficients of the expansions.  Smoothness is controlled
%  by controlling the size of an integrated squared linear differential
%  operator.  The integral is multiplied by a smoothing or bandwidth
%  parameter.
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
%  WTVEC    ... A vector of N weights, set to one by default, that can
%               be used to differentially weight observations.
%               WTVEC may also be an order N matrix (usually symmetric)
%               that also allows for correlations among observations.
%  FDNAMES  ... A cell array of length 3 with names for
%               1. argument domain, such as 'Time'
%               2. replications or cases
%               3. the function.
%  DFFACTOR ... A multiplier of df in GCV, set to one by default
%
%  Returned objects are:
%
%  FDOBJ   ... an object of class fd containing coefficients.
%  DF      ... a degrees of freedom measure.
%  GCV     ... a measure of lack of fit discounted for df.
%                 If the function is univariate, GCV is a vector 
%                 containing the error  sum of squares for each 
%                 function, and if the function is multivariate, 
%                 GCV is a NVAR by NCURVES matrix.
%  SSE     ... the error sums of squares.  
%                 SSE is a vector or matrix of the same size as 
%                 GCV.
%  PENMAT  ... the penalty matrix, if computed, otherwise [].
%  Y2CMAP  ... the matrix mapping the data to the coefficients.
%  ARGVALS ... the input set of argument values.
%  Y       ... the input array containing values of curves

%  Last modified 20 July 2010 by Jim Ramsay

if nargin < 3
    error('There is not at least three arguments.');
end

%  check ARGVALS

[argvals, n] = argcheck(argvals);

%  check Y

[y, ncurve, nvar, ndim] = ycheck(y, n);

%  check FDPAROBJ and get FDOBJ and LAMBDA

fdParobj = fdParcheck(fdParobj);
fdobj    = getfd(fdParobj);
lambda   = getlambda(fdParobj);
Lfdobj   = getLfd(fdParobj);

%  check LAMBDA

if lambda < 0
    lambda = 0;
end

%  set default argument values

if nargin < 6, dffactor = 1;        end
if nargin < 5
    fdnames = cell(1,3);
    fdnames{1} = 'arguments';
    fdnames{2} = 'replications';
    fdnames{3} = 'variables';
end

if nargin < 4, wtvec  = ones(n,1);  end

%  get BASIS and NBASIS

basisobj = getbasis(fdobj);
nbasis   = getnbasis(basisobj) - length(getdropind(basisobj));

%  check WTVEC

[wtvec, onewt, matwt] = wtcheck(n, wtvec);

%  check DFFACTOR

if isempty(dffactor), dffactor = 1.0;  end

%  ------------------------------------------------------------------
%                set up the linear equations for smoothing
%  ------------------------------------------------------------------

%  set up matrix of basis function values

basismat  = eval_basis(argvals, basisobj);

if n >= nbasis || lambda > 0
    
    %  The following code is for the coefficients completely determined
    
    if matwt
        basisw = basismat*wtvec;
    else
        basisw = basismat .* (wtvec * ones(1,nbasis));
    end
    Bmat   = basisw' * basismat;
    Bmat0  = Bmat;
    
    %  set up right side of equation
    
    if ndim < 3
        Dmat = basisw' * y;
    else
        Dmat = zeros(nbasis,ncurve,nvar);
        for ivar = 1:nvar
            Dmat(:,:,ivar) = basisw' * y(:,:,ivar);
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
        Bmat = Bmat + lambda .* penmat;
    else
        penmat = [];
    end
    
    %  compute inverse of Bmat
    
    if is_diag(Bmat)
        Bmatinv = diag(1./diag(Bmat));
    else
        Bmatinv = inv(Bmat);
    end
    
    %  ------------------------------------------------------------------
    %       Compute the coefficients defining the smooth and 
    %            summary properties of the smooth
    %  ------------------------------------------------------------------

    %  compute map from y to c
    
    y2cMap = Bmatinv * basisw';
    
    %  compute degrees of freedom of smooth
    
    df = full(sum(diag(Bmatinv * Bmat0)));
    
    %  solve normal equations for each observation
    
    if ndim < 3
        coef = Bmatinv * Dmat;
    else
        coef = zeros(nbasis, ncurve, nvar);
        for ivar = 1:nvar
            coef(:,:,ivar) = Bmatinv * Dmat(:,:,ivar);
        end
    end
    
else
    error(['The number of basis functions exceeds the number of ', ...
           'points to be smoothed.']);
    
end

%  ------------------------------------------------------------------
%            compute SSE, yhat, GCV and other fit summaries
%  ------------------------------------------------------------------

%  compute error sum of squares

if ndim < 3
    yhat = basismat * coef;
    SSE = sum((y - yhat).^2);
else
    SSE = zeros(nvar,ncurve);
    for ivar = 1:nvar
        coefi = squeeze(coef(:,:,ivar));
        yhati = basismat * coefi;
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

fdobj = fd(coef, basisobj, fdnames);


