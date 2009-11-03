function [fdobj, beta, df, gcv, SSE, penmat, y2cMap, argvals, y] = ...
    smooth_basis_LS(argvals, y, fdParobj, varargin)
%SMOOTH_BASIS_LS  Smooths data by penalized least squares.  The smooth 
%  curves are expressed as a basis function expansions, and this function 
%  computes the coefficients of the expansions.  Smoothness is controlled
%  by controlling the size of an integrated squared linear differential
%  operator.  The integral is multiplied by a smoothing or bandwidth
%  parameter.
%
%  In addition, an optional covariate or design matrix can be supplied,
%  having a row corresponding to each value in ARGVALS and Q columns,
%  Q being the number of covariates.  See the optional COVARIATES
%  parameter described below.
%
%  This version of function smooth_basis sets up the smoothing problem
%  as a least squares fitting problem, with the penalty term set up as
%  the smooth of a basis system toward 0.
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
%   A function call of the form
%             smooth_basis_LS(argvals, y, fdParobj,
%                             'PARAM1',val1,'PARAM2',val2,...)
%   can be used to input a number of optional parameters.  The first
%   of each pair of arguments is the name of the parameter, supplied
%   as a string in quotes, and the second is the value of the argument,
%   supplied as an object of the class required by the parameter.
%  
%   These optional parameters are:
%
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
%
%  Returned objects are:
%
%  FDOBJ   ... an object of class fd containing coefficients.
%  BETA    ... the regression coefficients for the covariates if supplied
%              or empty otherwise
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

%  Last modified 25 July 2010 by Jim Ramsay

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

%  set up default fdnames

deffdnames = cell(1,3);
deffdnames{1} = 'arguments';
deffdnames{2} = 'replications';
deffdnames{3} = 'variables';

% Process optional name/value pairs.

%  set up the argument names and argument default values

paramNames = {     'weight' 'fdnames' 'covariates' 'dfscale'};
paramDflts = {[]     deffdnames      []        1.0};

%  Use function getargs to obtain the values that are supplied, and
%  return the default values for those that are not.

[errid,errmsg,wtvec,fdnames,covariates,dfscale] = ...
         internal.stats.getargs(paramNames, paramDflts, varargin{:});
if ~isempty(errid)
    error(sprintf('smooth_basis_LS:optionargs:%s',errid),errmsg);
end

%  get BASIS and NBASIS

basisobj = getbasis(fdobj);
nbasis   = getnbasis(basisobj) - length(getdropind(basisobj));

%  check WTVEC

[wtvec, onewt] = wtcheck(n, wtvec);

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
    
    %  Multiply the basis matrix and the data pointwise by the square root 
    %  of the weight vector if the weight vector is not all ones.
    
    if ~onewt
        rtwtvec  = sqrt(wtvec);
        basismat = basismat .* repmat(rtwtvec,1,nbasis);
        rtwtmat  = repmat(rtwtvec,1,ncurve);
        if ndim < 3
            y = y .* rtwtmat;
        else
            for ivar=1:nvar
                y(:,:,ivar) = y(:,:,ivar).rtwtmat;
            end
        end
    end
    
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
        if nderiv > 0 && log10(Dsort(eiglow+1)/Dsort(1)) > -1e-12
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
            y = [y; zeros(nbasis-nderiv,ncurve)];
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
        if ~onewt
            basismat0(ind1,ind2) = covariates.*repmat(rtwtvec,1,q);
            basismat(ind1,ind2)  = covariates.*repmat(rtwtvec,1,q);
        else
            basismat0(ind1,ind2) = covariates;
            basismat(ind1,ind2)  = covariates;
        end
        if sparsewrd
            basismat0 = sparse(basismat0);
            basismat  = sparse(basismat);
        end
    end
    
    %  solve the least squares problem using the QR decomposition with
    %  one iteration to improve accuracy
    
    [Q,R]= qr(basismat,0);
    
    if ndim < 3
        coef = R\(R'\(basismat'*y));
        res  = y - basismat*coef;
        err  = R\(R'\(basismat'*res));
        coef = coef + err;
    else
        coef = zeros(nbasis, ncurve, nvar);
        for ivar = 1:nvar
            yi = squeeze(y(:,:,ivar));
            coefi = R\(R'\(basismat'*yi));
            resi  = yi - basismat*coefi;
            erri  = R\(R'\(basismat'*resi));
            coefi = coefi + erri;
            coef(:,:,ivar) = coefi;
        end
    end
    
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


