function [fdobj, beta, df, gcv, SSE, dev, stats, ...
          penmat, y2cMap, argvals, y] = ...
                       smooth_basis_GLM(argvals, y, fdParobj, varargin)
%SMOOTH_GLM  Smooths discrete curve represented by basis function
%  expansions fit by penalized least squares.
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
%  Optional arguments are input in pairs:  the first element of the pair
%     is a string specifying the property that the argument value defines,
%     and the second element is the value of the argument
%
%     Valid property/value pairs include:
%
%     Property        Value
%     ----------------------------------------------------------------
%     weight          vector of the same length as the data vector to be
%                     smoothed, containing nonnegative weights to be 
%                     applied to the data values
%     fdnames         A cell array of length 3 with names for
%                       1. argument domain, such as 'Time'
%                       2. replications or cases
%                       3. the function.
%     covariates      A N by Q matrix Z of covariate values used to augment
%                     the smoothing function, where N is the number of
%                     data values to be smoothed and Q is the number of
%                     covariates.  The process of augmenting a smoothing 
%                     function in this way is often called "semi-parametric 
%                     regression".  The default is the empty object [].
%     dfscale         A scalar value multiplying the degrees of freedom
%                     in the definition of the generalized 
%                     cross-validated or GCV criterion for selecting the
%                     bandwidth parameter LAMBDA.  It was recommended by
%                     Chong Gu that this be a number slightly larger than
%                     1.0, such as 1.2, to prevent under-smoothing,
%                     The default is 1.0.
%     family          a character string containing one of:
%                       'normal'  
%                       'binomial'
%                       'poisson'
%                       'gamma'
%                       'inverse gaussian'
%                     the value determines which of the link functions in
%                     the generalized linear model (GLM) family is to be
%                     used.  The default is 'normal'.
%      control        a struct object controlling iterations with members:
%                       epsilon  convergence criterion (default 1e-8)
%                       maxit    max. iterations       (default 25)
%                       trace    output iteration info (0)
%      start          a vector containing starting values for coefficients
%                      
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

%  Last modified 26 July 2010 by Jim Ramsay

if nargin < 3
    error('There is not at least three arguments.');
end

%  check ARGVALS

[argvals, n] = argcheck(argvals);

%  check Y

% [y, ncurve, nvar, ndim] = ycheck(y, n);
ydim = size(y);
if length(ydim) == 2 && ydim(1) == 1
    y = y';
end

ydim = size(y);  %  number of observations
if ydim(1) ~= n
    error('Y is not the same length as ARGVALS.');
end

%  set number of curves and number of variables

sizey = size(y);
ndim  = length(sizey);
switch ndim
    case 1
        ncurve = 1;
        nvar   = 1;
    case 2
        ncurve = sizey(2);
        nvar   = 1;
    case 3
        ncurve = sizey(2);
        nvar   = sizey(3);
    otherwise
        error('Second argument must not have more than 3 dimensions');
end

%  check FDPAROBJ and get FDOBJ and LAMBDA

fdParobj = fdParcheck(fdParobj);
fdobj    = getfd(fdParobj);
lambda   = getlambda(fdParobj);
Lfdobj   = getLfd(fdParobj);

%  check LAMBDA

if lambda < 0, lambda = 0;  end

%  get BASIS and NBASIS

basisobj = getbasis(fdobj);
nbasis   = getnbasis(basisobj) - length(getdropind(basisobj));

%  set default argument values

deffdnames = cell(1,3);
deffdnames{1} = 'arguments';
deffdnames{2} = 'replications';
deffdnames{3} = 'variables';

% Process optional name/value pairs.

%  set up the argument names and argument default values

paramNames = {     'weight' 'fdnames' 'covariates' 'dfscale' 'family'};
paramDflts = {[]     deffdnames      []        1.0     'normal'};

%  Use function getargs to obtain the values that are supplied, and
%  return the default values for those that are not.

[errid,errmsg,wtvec,fdnames,covariates,dfscale,family] = ...
         internal.stats.getargs(paramNames, paramDflts, varargin{:});
if ~isempty(errid)
    error(sprintf('smooth_basis_LS:optionargs:%s',errid),errmsg);
end

%  check WTVEC

[wtvec, onewt] = wtcheck(n, wtvec);
if onewt
    wtvec = ones(n,1);
end

%  check FDNAMES

if ~iscell(fdnames)
    error('smooth_basis_LS:fdnames', ...
          'Optional argument FDNAMES is not a cell array.');
end

if length(fdnames) ~= 3
    error('smooth_basis_LS:fdnames', ...
          'Optional argument FDNAMES is not of length 3.');
end

%  check COVARIATES

q = 0;
if ~isempty(covariates)
    if ~isnumeric(covariates)
        error('smooth_basis_LS:covariates', ...
            'Optional argument COVARIATES is not numeric.');
    end
    if size(covariates,1) ~= n
        error('smooth_basis_LS:covariates', ...
            'Optional argument COVARIATES has incorrect number of rows.');
    end
    q = size(covariates,2);
end

%  ------------------------------------------------------------------
%                set up the linear equations for smoothing
%  ------------------------------------------------------------------

%  set up matrix of basis function values

basismat  = eval_basis(argvals, basisobj);

if n >= nbasis || lambda > 0
    
    %  The following code is for the coefficients completely determined
    
    %  set up additional rows of the least squares problem for the
    %  penalty term.

    basismat0 = basismat;
    y0        = y;
    
    if lambda > 0
        nderiv  = getnderiv(Lfdobj);
        penmat  = eval_penalty(basisobj, Lfdobj);
        [V,D] = eig(full(penmat));
        Dvec  = diag(D);
        [Dsort, Isort] = sort(Dvec, 'descend');
        Vsort = V(:,Isort);
        %  Check that the lowest eigenvalue in the series that is to be
        %  kept is positive.
        eiglow = nbasis - nderiv;
        if Dsort(eiglow) <= 0
            error('smooth_basis:eig', ...
                  ['Eigenvalue(NBASIS-NDERIV) of penalty matrix ', ...
                   'is not positive; check penalty matrix.']);
        end
        %  Check that the highest eigenvalue that is not used is small
        %  relative to the largest eigenvalue.
       
        if nderiv > 0 && Dsort(eiglow+1) > 0 && ...
                         log10(Dsort(eiglow+1)/Dsort(1)) > -1e-12
            error('smooth_basis:eig', ...
                  ['Eigenvalue(NBASIS-NDERIV+1) of penalty matrix ', ...
                   'is not small relative to eigenvalue(1); ', ...
                   'check penalty matrix.']);
        end
        %  Compute the square root of the penalty matrix in the subspace
        %  spanned by the first N - NDERIV eigenvectors
        ind = 1:eiglow;
        penfac = Vsort(:,ind)*diag(sqrt(Dsort(ind)));
        %  Augment basismat by sqrt(lambda).*penfac'
        basismat = [basismat; sqrt(lambda).*penfac'];
        %  Augment data vector by n - nderiv 0's
        if ndim < 3
            y     = [y; zeros(nbasis-nderiv,ncurve)];
            wtvec = [wtvec; ones(nbasis-nderiv,1)];
        else
            for ivar=1:nvar
                y(:,:,ivar) = [y(:,:,ivar); zeros(n-nderiv,ncurve)];
            end
        end
    end
    
    %  augment BASISMAT0 and BASISMAT by the covariate matrix 
    %  if it is supplied
    
    if ~isempty(covariates)
        ind1 = 1:n;
        ind2 = (nbasis+1):(nbasis+q);
        sparsewrd = issparse(basismat0);
        basismat0 = full([basismat0, zeros(size(basismat0,1),q)]);
        basismat  = full([basismat,  zeros(size(basismat, 1),q)]);
        basismat0(ind1,ind2) = covariates;
        basismat(ind1,ind2)  = covariates;
        if sparsewrd
            basismat0 = sparse(basismat0);
            basismat  = sparse(basismat);
        end
    end
    
    %  define data indicator array DATAIND
    
    dataind = zeros(size(basismat,1),1);
    dataind(1:n) = 1;
    
    %  ------------------------------------------------------------------
    %               compute solution using function GLMFIT_FDA
    %  ------------------------------------------------------------------

%     [bb,dev,stats] = glmfit_fda(x,y,distr,varargin)
%  set up for going into GLMFIT_FDA
    x = basismat;
    distr = family;
    varargin = cell(6,1);
    varargin{1} = 'constant';
    varargin{2} = 'off';
    varargin{3} = 'weights';
    varargin{4} = wtvec;
    varargin{5} = 'dataind';
    varargin{6} = dataind;
    
    if ndim < 3
        coef  = zeros(nbasis+q,ncurve);
        dev   = zeros(ncurve,1);
        stats = cell(ncurve,1);
        for icurve=1:ncurve
            [coefi,devi,statsi] = glmfit_fda(basismat, y, family, ...
                  'constant', 'off', 'weights', wtvec, 'dataind', dataind);
            coef(:,icurve) = coefi;
            dev(icurve)    = devi;
            stats{icurve}  = statsi;
        end
    else
        coef  = zeros(nbasis+1,ncurve,nvar);
        dev   = zeros(ncurve,nvar);
        stats = cell(ncurve,nvar);
        for ivar=1:nvar
            for icurve=1:ncurve
                [coefi,devi,statsi] = glmfit_fda(basismat, y, family, ...
                  'constant', 'off', 'weights', wtvec, 'dataind', dataind);
                coef(:,icurve,ivar) = coefi;
                dev(icurve,ivar)    = devi;
                stats{icurve,ivar}  = statsi;
            end
        end
    end
    
    %  Compute Q-R decomposition of extended basis matrix
            
    [Q,R]= qr(basismat,0);
    
    %  compute basismat*R^{-1}
    
    MapFac = R'\basismat0';
    
    %  compute map from y to c
    
    y2cMap = R\MapFac;
    
    %  compute degrees of freedom of smooth
    
    df = full(sum(diag(MapFac*MapFac')));
 
else
    error(['The number of basis functions exceeds the number of ', ...
           'points to be smoothed.']);    
end

%  ------------------------------------------------------------------
%            compute SSE, yhat, GCV and other fit summaries
%  ------------------------------------------------------------------

%  compute error sum of squares

if ndim < 3
    yhat = basismat0 * coef;
    SSE  = sum((y0 - yhat).^2);
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
    gcv = (SSE./n)./((n - dfscale*df)/n)^2;
else
    gcv = NaN;
end

%  set up the functional data object

if ndim < 3
    fdobj = fd(coef(1:nbasis,:),   basisobj, fdnames);
else
    fdobj = fd(coef(1:nbasis,:,:), basisobj, fdnames);
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






