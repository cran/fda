function [fdobj, df, gcv, beta, SSE, penmat, y2cMap, argvals, y] = ...
                      smooth_basis(argvals, y, fdParobj, varargin)
%SMOOTH_BASIS  Smooths data by penalized least squares.  The smooth 
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
%  ARGVALS  ... A set of n argument values, set by default to be 
%               equally spaced on the unit interval (0,1).
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
%     weight          vector of length n, containing nonnegative weights  
%                     to be applied to the data values, or
%                     a symmetric positive definite matrix of order n
%                     also w, wt, wgt, weights
%                     See note below about a change in the roughness 
%                     penalty related to using weights
%     fdnames         A cell array of length 3 with names for
%                       1. argument domain, such as 'Time'
%                       2. replications or cases
%                       3. the function.
%                     also f, fdname
%     covariates      A N by Q matrix Z of covariate values used to augment
%                     the smoothing function, where N is the number of
%                     data values to be smoothed and Q is the number of
%                     covariates.  The process of augmenting a smoothing 
%                     function in this way is often called "semi-parametric 
%                     regression".  The default is the empty object [].
%                     also c, cov, covariate
%     method          The method for computing coefficients.  The usual
%                     method computes cross-product matrices of the basis
%                     value matrix, adds the roughness penalty, and uses
%                     the Choleski decomposition of this to compute
%                     coefficients, analogous to using the normal equations
%                     in least squares fitting.  But this approach, while
%                     fast, contributes unnecessary rounding error, and the
%                     qr decomposition of the augmented basis matrix is
%                     prefererable.   But nothing comes for free, and the
%                     computational overhead of the qr approach can be a
%                     serious problem for large problems (n of 1000 or
%                     more).  For this reason, the default is 'method' ==
%                     'chol', but if 'method' == 'qr', the qr decomposition
%                     is used.
%     dfscale         A scalar value multiplying the degrees of freedom
%                     in the definition of the generalized 
%                     cross-validated or GCV criterion for selecting the
%                     bandwidth parameter LAMBDA.  It was recommended by
%                     Chong Gu that this be a number slightly larger than
%                     1.0, such as 1.2, to prevent under-smoothing,
%                     The default is 1.0.
%                     also d, df
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
%  BETA    ... the regression coefficients for the covariates if supplied
%              or empty otherwise
%  SSE     ... the error sums of squares.  
%                 SSE is a vector or matrix of the same size as 
%                 GCV.
%  PENMAT  ... the penalty matrix, if computed, otherwise [].
%  Y2CMAP  ... the matrix mapping the data to the coefficients.
%  ARGVALS ... the input set of argument values.
%  Y       ... the input array containing values of curves

%  Last modified 23 March 2012 by Jim Ramsay

%  Note: (21 December 2011) The factor sum(wtvec)/n for matwt == 0
%  or trace(wtvec)/n if matwt ~= 0 has been added to the multiplier
%  of the roughness penalty so as to make the total penalized least
%  squares homogeneous with respect to a scalar multiplier of the
%  weights.  That is, if wtvec is multiplied by constant C, then so
%  is the roughness penalty, and the PLS therefore by C as well.
%  This important because otherwise the amount of smoothing would
%  vary if the weights were changed by a scale factor.

%  ---------------------------------------------------------------------
%                      Check argments
%  ---------------------------------------------------------------------

if nargin < 3
    error('There is not at least three arguments.');
end

%  check ARGVALS

[argvals, n] = argcheck(argvals);

%  check Y

[y, ncurve, nvar, ndim] = ycheck(y, n);
y0 = y;  %  preserve a copy of y;

%  check FDPAROBJ and get FDOBJ and LAMBDA

fdParobj = fdParcheck(fdParobj);
fdobj    = getfd(fdParobj);
lambda   = getlambda(fdParobj);
Lfdobj   = getLfd(fdParobj);
penmat   = getpenmat(fdParobj);

%  check LAMBDA

if lambda < 0
    lambda = 0;
end

%  set up default fdnames

deffdnames = cell(1,3);
deffdnames{1} = 'arguments';
deffdnames{2} = 'replications';
deffdnames{3} = 'variables';

% Which style of calling sequence is being used:  
%    name -- value pairs or fixed order?

NameValStyle = true;
if ~isempty(varargin) && nargin <= 7
   va1 = varargin{1};
   if ~ischar(va1) || isempty(va1)
      NameValStyle = false;
   end
end

%  set default argument values

wtvec      = [];
fdnames    = deffdnames;
covariates = [];
method     = 'chol';
dfscale    = 1;

if NameValStyle
    
    % Process optional number of name -- value pairs.
    
    nargpr = nargin - 3;
    if floor(nargpr/2)*2 ~= nargpr
        error(['The number of argments after the first three ', ...
            'is not an even number.']);
    end
    
    for ipr=5:2:nargin
        ArgName = varargin{ipr-4};
        if     strcmp(ArgName, 'w')         || ...
               strcmp(ArgName, 'wt')        || ...
               strcmp(ArgName, 'wgt')       || ...
               strcmp(ArgName, 'weight')    || ...
               strcmp(ArgName, 'weights')
            wtvec      = varargin{ipr-3};
        elseif strcmp(ArgName, 'f')         || ...
               strcmp(ArgName, 'fdname')    || ...
               strcmp(ArgName, 'fdnames')
            fdnames    = varargin{ipr-3};
        elseif strcmp(ArgName, 'c')         || ...
               strcmp(ArgName, 'cov')       || ...
               strcmp(ArgName, 'covariate') || ...
               strcmp(ArgName, 'covariates')
            covariates = varargin{ipr-3};
        elseif strcmp(ArgName, 'm')         || ...
               strcmp(ArgName, 'meth')      || ...
               strcmp(ArgName, 'method')
            method = varargin{ipr-3};
        elseif strcmp(ArgName, 'd')         || ...
               strcmp(ArgName, 'df')        || ...
               strcmp(ArgName, 'dfscl')     || ...
               strcmp(ArgName, 'dfscale')
            dfscale    = varargin{ipr-3};
        else
            error('An argument name is unrecognizable.');
        end
    end
else
    
    %  process argument values in fixed order
    
    if nargin >= 4,  wtvec      = varargin{1};   end
    if nargin >= 5,  fdnames    = varargin{2};   end
    if nargin >= 6,  covariates = varargin{3};   end
    if nargin >= 7,  method     = varargin{4};   end
    if nargin == 8,  dfscale    = varargin{5};   end
    if nargin > 7
        error('More than seven non-named arguments found.');
    end
end

%  get BASIS and NBASIS

basisobj = getbasis(fdobj);
nbasis   = getnbasis(basisobj) - length(getdropind(basisobj));
ind      = 1:nbasis;

%  check WTVEC

[wtvec, onewt, matwt] = wtcheck(n, wtvec);

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

if strcmp(method, 'chol')
    
    %  -----------------------------------------------------------------
    %  use the default choleski decomposition of the crossproduct of the
    %  basis value matrix plus the roughness penalty
    %  -----------------------------------------------------------------
    
    if n >= nbasis || lambda > 0
        
        %  augment BASISMAT0 and BASISMAT by the covariate matrix
        %  if it is supplied
        
        if ~isempty(covariates)
            ind1 = 1:n;
            ind2 = (nbasis+1):(nbasis+q);
            sparsewrd = issparse(basismat);
            basismat  = full([basismat,  zeros(size(basismat, 1),q)]);
            basismat(ind1,ind2)  = covariates;
            if sparsewrd
                basismat  = sparse(basismat);
            end
        end
        
        %  The following code is for the coefficients completely determined
        
        if matwt
            wtmat = wtvec;
            wtfac = chol(wtmat);
            basisw = wtmat*basismat;
        else
            rtwtvec  = sqrt(wtvec);
            rtwtmat  = repmat(rtwtvec,1,ncurve);
            basisw   = basismat .* (wtvec * ones(1,nbasis+q));
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
            if isempty(penmat)
                penmat = eval_penalty(basisobj, Lfdobj);
            end
            Bnorm   = sqrt(sum(sum(Bmat.^2)));
            pennorm = sqrt(sum(sum(penmat.^2)));
            condno  = pennorm/Bnorm;
            if lambda*condno > 1e12
                lambda = 1e12/condno;
                warning('Wid2:reduce', ...
                    ['LAMBDA reduced to ',num2str(lambda), ...
                    ' to prevent overflow']);
            end
            if ~isempty(covariates)
                penmat = [penmat,          zeros(nbasis,q); ...
                          zeros(q,nbasis), zeros(q,q)];
            end
            if onewt
                smoothfac = lambda;
            else
                if matwt
                    smoothfac = trace(wtvec)*lambda/n;
                else
                    smoothfac =   sum(wtvec)*lambda/n;
                end
            end
            Bmat = Bmat + smoothfac .* penmat;
        else
            penmat = zeros(nbasis);
            Bmat   = Bmat0;
        end
        
        %  compute Choleski factor of Bmat
        
        [Bmatfac, rnkp1] = chol(Bmat);
        if rnkp1 ~= 0
            error(['The matrix BMAT  + LAMBDA .* PENMAT ', ...
                'is not positive definite.']);
        end
        
        %  compute coefficient matrix or array
        
        if ndim < 3
            coef = Bmatfac\(Bmatfac'\Dmat);
        else
            coef = zeros(nbasis, ncurve, nvar);
            for ivar = 1:nvar
                coef(:,:,ivar) = Bmatfac\(Bmatfac'\Dmat(:,:,ivar));
            end
        end
        
    else
        error(['The number of basis functions exceeds the number of ', ...
            'points to be smoothed.']);
        
    end
    
elseif strcmp(method, 'qr')
    
    %  -------------------------------------------------------------
    %  computation of coefficients using the qr decomposition of the
    %  augmented basis value matrix
    %  -------------------------------------------------------------
    
    if n >= nbasis || lambda > 0
        
        %  The following code is for the coefficients completely determined
        
        %  Multiply the basis matrix and the data pointwise by 
        %  the square root of the weight vector if the weight vector 
        %  is not all ones.  If the weights are in a matrix, 
        %  multiply the basis matrix by its Choleski factor.
        
        if ~onewt
            if matwt
                wtmat = (wtvec + t(wtvec))/2;
                wtfac = chol(wtvmat);
                basismat = wtfac * basismat;
                if ndim < 3
                    y = wtfac * y;
                else
                    for ivar=1:nvar
                        y(:,:,ivar) = wtfac * squeeze(y(:,:,ivar));
                    end
                end
            else
                rtwtvec  = sqrt(wtvec);
                basismat = basismat .* repmat(rtwtvec,1,nbasis);
                rtwtmat  = repmat(rtwtvec,1,ncurve);
                if ndim < 3
                    y = rtwtmat .* y;
                else
                    for ivar=1:nvar
                        y(:,:,ivar) = rtwtmat .* y(:,:,ivar);
                    end
                end
            end
        end
        
        %  set up additional rows of the least squares problem for the
        %  penalty term.
        
        if lambda > 0
            nderiv = getnderiv(Lfdobj);
            if isempty(penmat)
                penmat = eval_penalty(basisobj, Lfdobj);
            end
            [V,D]  = eig(full(penmat));
            Dvec   = diag(D);
            [Dsort, Isort] = sort(Dvec, 'descend');
            Vsort  = V(:,Isort);
            %  Check that the lowest eigenvalue in the series that is to be
            %  kept is positive.
            eiglow = nbasis - nderiv;
            if Dsort(eiglow) <= 0
                error('smooth_basis:eig', ...
                    ['Eigenvalue(NBASIS-NDERIV) of penalty matrix ', ...
                    'is not positive; check penalty matrix.']);
            end
            %  Compute the square root of the penalty matrix in the subspace
            %  spanned by the first N - NDERIV eigenvectors
            indeig = 1:eiglow;
            penfac = Vsort(:,indeig)*diag(sqrt(Dsort(indeig)));
            %  Augment basismat by sqrt(lambda).*penfac'
            if onewt
                smoothfac = sqrt(lambda);
            else
                if matwt
                    smoothfac = sqrt(trace(wtvec)*lambda/n);
                else
                    smoothfac = sqrt(  sum(wtvec)*lambda/n);
                end                    
            end
            basismat_aug = [basismat; smoothfac.*penfac'];
            %  Augment data vector by n - nderiv 0's
            if ndim < 3
                y = [y; zeros(nbasis-nderiv,ncurve)];
            else
                for ivar=1:nvar
                    y(:,:,ivar) = [y(:,:,ivar); zeros(n-nderiv,ncurve)];
                end
            end
        else
            penmat = [];
        end
        
        %  augment BASISMAT0 and BASISMAT by the covariate matrix
        %  if it is supplied
        
        if ~isempty(covariates)
            ind1 = 1:n;
            ind2 = (nbasis+1):(nbasis+q);
            sparsewrd = issparse(basismat);
            basismat_aug = ...
                full([basismat_aug, zeros(size(basismat, 1),q)]);
            if ~onewt
                if matwt
                    basismat_aug(ind1,ind2)  = ...
                        wtfac*covariates;
                else
                    basismat_aug(ind1,ind2)  = ...
                        covariates.*repmat(rtwtvec,1,q);
                end
            else
                basismat(ind1,ind2)  = covariates;
            end
            if sparsewrd
                basismat_aug = sparse(basismat_aug);
            end
            if ~isempty(covariates)
                penmat = [penmat,          zeros(nbasis,q); ...
                          zeros(q,nbasis), zeros(q,q)];
            end
        end
        
        %  solve the least squares problem using the QR decomposition
        
        [~,R]= qr(basismat_aug,0);
        
        if ndim < 3
            coef = R\(R'\(basismat'*y));
            res  = y - basismat*coef;
            err  = R\(R'\(basismat'*res));
            coef = coef + err;
            if ~isempty(covariates)
                beta = coef(ind2,:);
                coef = coef(ind1,:);
            end
        else
            coef = zeros(nbasis, ncurve, nvar);
            beta = zeros(q, ncurve, nvar);
            for ivar = 1:nvar
                yi = squeeze(y(:,:,ivar));
                coefi = R\(R'\(basismat'*yi));
                resi  = yi - basismat*coefi;
                erri  = R\(R'\(basismat'*resi));
                coefi = coefi + erri;
                if ~isempty(covariates)
                    betai = coefi(ind2,:);
                    coefi = coefi(ind1,:);
                end
                coef(:,:,ivar) = coefi;
                beta(:,:,ivar) = betai;
            end
        end
        
    else
        error(['The number of basis functions exceeds the number of ', ...
            'points to be smoothed and lambda is 0.']);
    end

else
    error('Unrecognizable value of method.');
end

if isempty(covariates)
    fdobj = fd(coef, basisobj, fdnames);
else
    if ndim < 3
        fdobj = fd(coef(ind,:),   basisobj, fdnames);
    else
        fdobj = fd(coef(ind,:,:), basisobj, fdnames);
    end        
end

if nargout == 1
    return;
end

%  ------------------------------------------------------------------
%            compute SSE, yhat, GCV and other fit summaries
%  ------------------------------------------------------------------

%  compute map from y to c

if onewt
    temp   = basismat'*basismat;
    if ~isempty(penmat)
        temp = temp + lambda.*penmat;
    end
    L      = chol(temp);
    MapFac = L'\basismat';
    y2cMap = L\MapFac;
else
    if ~matwt
        wtmat  = sparse(diag(wtvec));
    end
    temp   = basismat'*wtmat*basismat;
    if ~isempty(penmat)
        temp = temp + lambda.*penmat;
    end
    L      = chol(temp);
    MapFac = L'\basismat';
    y2cMap = L\MapFac*wtmat;
end

%  compute degrees of freedom of smooth

df = trace(y2cMap*basismat);

%  compute error sum of squares

if ndim < 3
    yhat = basismat * coef;
    if onewt
        SSE = sum((y0 - yhat).^2);
    else
        if matwt
            SSE = sum((wtfac*(y0 - yhat)).^2);
        else
            SSE = sum((rtwtmat.*(y0 - yhat)).^2);
        end
    end
else
    SSE = zeros(nvar,ncurve);
    for ivar = 1:nvar
        coefi = squeeze(coef(:,:,ivar));
        yhati = basismat * coefi;
        yi    = squeeze(y(1:n,:,ivar));
        SSE(ivar,:) = sum((yi - yhati).^2);
    end
end

%  compute  GCV index

if df < n
    gcv = (SSE./n)./((n - dfscale*df)/n)^2;
else
    gcv = NaN;
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

%  restore penmat to its original if not empty and covariates used

if ~isempty(penmat) && ~isempty(covariates)
    penmat = penmat(1:nbasis,1:nbasis);
end
