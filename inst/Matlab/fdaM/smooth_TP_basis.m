function [TP_fdobj, df, gcv, coefmat, SSE, y2cMap] = ...
    smooth_TP_basis(sarg, targ, y, TP_basisobj, ...
                    sLfd, tLfd, slambda, tlambda, fdnames)
%SMOOTH_TP_BASIS  Smooths discrete surface values over a rectangular 
%  lattice to estimate a smoothing function f(s,t) 
%   using penalized basis functions.
%  Arguments for this function:
%
%  sarg     ... A set of argument values for the row dimension s. 
%  targ     ... A set of argument values for the col dimension t. 
%  Y        ... an array containing surface values
%  FDPARS   ... A functional parameter or fdPar object for 
%               variation along the row dimension s.
%  FDPART   ... A functional parameter or fdPar object for 
%               variation along the col dimension t.
%  FDNAMES  ... A cell of length 3 with names for
%               1. argument s
%               2. argument t
%               3. the function f(s,t).
%  Returns:
%    FDOBJ  ...  an object of class fd containing coefficients.
%    DF     ...  a degrees of freedom measure.
%    GCV    ...  a measure of lack of fit discounted for df.
%                If the function is univariate, GCV is a vector 
%                containing the error  sum of squares for each 
%                function, and if the function is multivariate, 
%                GCV is a NVAR by NCURVES matrix.
%    COEF   ...  the coefficient matrix for the basis function
%                  expansion of the smoothing function
%    SSE    ...  the error sums of squares.  
%                SSE is a vector or matrix of the same size as 
%                GCV.
%    PENMAT ...  the penalty matrix.
%    Y2CMAP ...  the matrix mapping the data to the coefficients.

%  Last modified 9 March 2011

if nargin < 4
    error('There is not at least four arguments.');
end

%  set default argument values

if nargin < 9
    fdnames{1} = 'argument s';
    fdnames{2} = 'argument s';
    fdnames{3} = 'function';
end

if nargin < 8, tlambda = 0;  end
if nargin < 7, slambda = 0;  end

if nargin < 6, tLfd = [];    end
if nargin < 5, sLfd = [];    end

%  check argument values

[sarg, ns] = argcheck(sarg);
[targ, nt] = argcheck(targ);

%  check Y

if ~strcmp(class(y), 'double')
    error('Y is not of class double.');
end

ydim = size(y);  %  number of observations
if length(ydim) ~= 2
    error('Y is not a matrix.');
end
if ydim(1) ~= ns
    error('Number of rows of Y is not the same length as SARG.');
end
if ydim(2) ~= nt
    error('Number of columns of Y is not the same length as TARG.');
end

yvec = reshape(y, ns*nt, 1);

%  check FDPARS, FDPART and BASES, LBFDOBJ'S and LAMBDA'S

TPbasisStruct = getbasispar(TP_basisobj);
sbasis = TPbasisStruct.basis1;
tbasis = TPbasisStruct.basis2;

snbasis = getnbasis(sbasis);
tnbasis = getnbasis(tbasis);

%  check LAMBDA's

if slambda < 0
    warning('Value of LAMBDAS was negative, 0 used instead.');
    slambda = 0;
end
if tlambda < 0
    warning('Value of LAMBDAT was negative, 0 used instead.');
    tlambda = 0;
end

%  ------------------------------------------------------------------
%                set up the linear equations for smoothing
%  ------------------------------------------------------------------

%  set up matrix of basis function values

sbasismat = eval_basis(sarg, sbasis);
tbasismat = eval_basis(targ, tbasis);
basismat  = kron(tbasismat,sbasismat);

if ns*nt >= snbasis*tnbasis || slambda > 0 || tlambda > 0
    
    %  The following code is for the coefficients completely determined
    
    Bmat   = basismat'*basismat;
    Bmat0  = Bmat;
    
    %  set up right side of equation
    
    Dmat = basismat' * yvec;
    
    %  set up regularized cross-product matrix BMAT

    if slambda > 0 || ~isempty(sLfd)
        penmats = eval_penalty(sbasis, sLfd);
        Bnorm   = norm(Bmat,    'fro');
        pennorm = norm(penmats, 'fro');
        condno  = pennorm/Bnorm;
        if slambda*condno > 1e12
            slambda = 1e12/condno;
            warning('Wid2:reduce', ...
                ['LAMBDAS reduced to ',num2str(slambda), ...
                    ' to prevent overflow']);
        end
        Bmat = Bmat + slambda .* kron(eye(tnbasis),penmats);
    end
    
    if tlambda > 0 || ~isempty(tLfd)
        penmatt = eval_penalty(tbasis, tLfd);
        Bnorm   = norm(Bmat,    'fro');
        pennorm = norm(penmatt, 'fro');
        condno  = pennorm/Bnorm;
        if tlambda*condno > 1e12
            tlambda = 1e12/condno;
            warning('Wid2:reduce', ...
                ['LAMBDAS reduced to ',num2str(tlambda), ...
                    ' to prevent overflow']);
        end
        Bmat = Bmat + tlambda .* kron(penmatt, eye(snbasis));
    end

    %  compute inverse of Bmat
    
    Bmatinv = inv(Bmat);
    
    %  ------------------------------------------------------------------
    %       Compute the coefficients defining the smooth and 
    %            summary properties of the smooth
    %  ------------------------------------------------------------------

    %  compute map from y to c
    
    y2cMap = Bmatinv * basismat';
    
    %  compute degrees of freedom of smooth
    
    df = full(sum(diag(Bmatinv * Bmat0)));
    
    %  solve normal equations for each observation
    
    opts.SYM = true;
    coef = linsolve(full(Bmat),full(Dmat),opts);
%     coefmat = reshape(coef, snbasis, tnbasis);
    
else
    error(['The number of basis functions exceeds the number of ', ...
           'points to be smoothed.']);
    
end

%  ------------------------------------------------------------------
%            compute SSE, yhat, GCV and other fit summaries
%  ------------------------------------------------------------------

%  compute error sum of squares

yhat = basismat * coef;
SSE = sum((yvec - yhat).^2);

%  compute  GCV index

if df < ns*nt
    gcv = (SSE./ns/nt)./((ns*nt - df)/ns/nt)^2;
else
    gcv = NaN;
end

TP_fdobj = fd(coef, TP_basisobj, fdnames);


