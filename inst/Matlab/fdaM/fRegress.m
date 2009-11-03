function fRegressStruct = ...
        fRegress(yfdPar, xfdcell, betacell, wt, y2cMap, SigmaE)
%  FREGRESS  Fits a functional linear model using multiple 
%  functional independent variables with the dependency being 
%  pointwise or concurrent.  
%  The case of a scalar independent variable is included by treating
%  it as a functional independent variable with a constant basis
%  and a unit coefficient.
%
%  Arguments:
%  YFDPAR   ... an object for the dependent variable, 
%               which may be: 
%                   a functional data object, 
%                   a functional parameter (fdPar) object, or 
%                   a vector
%  XFDCELL  ... a cell object of length p with each cell 
%               containing an object for an independent variable.
%               the object may be: 
%                   a functional data object or
%                   a vector
%  BETACELL ... a cell object of length p with each cell
%               containing a functional parameter object for 
%               the corresponding regression function.
%  WT       ... a vector of N nonnegative weights for observations,
%               where N is the number of observations. 
%  Y2CMAP   ... the matrix mapping from the vector of observed values
%               to the coefficients for the dependent variable.
%               This is output by function SMOOTH_BASIS.  If this is
%               supplied, confidence limits are computed, otherwise not.
%  SIGMAE   ... Estimate of the covariances among the residuals.  This
%               can only be estimated after a preliminary analysis
%               with FREGRESS.
%  
%  Returns:  
%  FREGRESSSTRUCT...  A struct object containing members with names:
%    yfdPar      ... first  argument of FREGRESS
%    xfdcell     ... second argument of FREGRESS
%    betacell    ... third  argument of FREGRESS
%    betahat     ... estimated regression functions 
%    yhat        ... if yfdPar is a functional parameter object, this is a
%                    functional data object containing fitted functions
%                    if yfdPar is a vector, this is a vector containing
%                    fitted values
%    Cmat        ... coefficient matrix for the linear system defining
%                    the regression coefficient basis vector
%    Dmat        ... right side vector for the linear system defining
%                    the regression coefficient basis vector
%    Cmatinv     ... inverse of the coefficient matrix, needed for
%                    function FREGRESS_STDERR that computes standard errors
%    wt          ... weights for observations
%    df          ... degrees of freedom for fit
%    y2cMap      ... input matrix mapping y observations into y
%                    coefficients
%    SigmaE      ... input covariance of residuals
%    betastderr  ... estimated standard error functional data objects for
%                    for regression coefficiens
%    bvar        ... covariance matrix for coefficients defining regression
%                    functions
%    c2bMap      ... matrix mapping y coefficients into b coefficients

%  Last modified 26 July 2010 by Jim Ramsay

if nargin < 3
    error('Less than three arguments supplied.');
end

%  set default values for last two arguments

if nargin < 6,  SigmaE = [];  end
if nargin < 5,  y2cMap = [];  end
if nargin < 4,  wt     = [];  end

[yfdPar, xfdcell, betacell, wt, rangeval] = ...
               fRegress_argcheck(yfdPar, xfdcell, betacell, wt);
          
p = length(xfdcell);
N = size(getcoef(xfdcell{1}),2);

wtconstant = (var(wt) == 0);

onesbasis = create_constant_basis(rangeval);
onesfd    = fd(1, onesbasis);

%  --------------------------------------------------------------------------
%  branch depending on whether the dependent variable is functional or scalar 
%  --------------------------------------------------------------------------

if isa_fdPar(yfdPar)
    
    %  ----------------------------------------------------------------
    %              YFDPAR is a functional parameter object
    %  ----------------------------------------------------------------
    
    %  extract dependent variable information
    
    yfdobj    = getfd(yfdPar);
    ycoef     = getcoef(yfdobj);
    ybasisobj = getbasis(yfdobj);
    rangeval  = getbasisrange(ybasisobj);
    ynbasis   = getnbasis(ybasisobj);
    
    if length(size(ycoef)) > 2
        error('YFDOBJ from YFDPAR is not univariate.');
    end
    
    %  -------- set up the linear equations for the solution  ----------
    
    %  compute the total number of coefficients to be estimated
    
    ncoef = 0;
    for j=1:p
        betafdParj = betacell{j};
        if getestimate(betafdParj)
            ncoefj  = size(getcoef(getfd(betafdParj)),1);
            ncoef   = ncoef + ncoefj;
        end
    end
    
    Cmat = zeros(ncoef,ncoef);
    Dmat = zeros(ncoef,1);
    
    %  loop through rows of CMAT
    
    mj2 = 0;
    for j=1:p
        betafdParj = betacell{j};
        if getestimate(betafdParj)
            betafdj    = getfd(betafdParj);
            betabasisj = getbasis(betafdj);
            ncoefj     = length(getcoef(betafdj));
            %  row indices of CMAT and DMAT to fill
            mj1    = mj2 + 1;
            mj2    = mj2 + ncoefj;
            indexj = mj1:mj2;
            %  compute right side of equation DMAT
            xfdj = xfdcell{j};
            if wtconstant
                xyfdj = xfdj.*yfdobj;
            else
                xyfdj = (xfdj.*wt).*yfdobj;
            end
            wtfdj = sum(xyfdj);
            Dmatj = inprod(betabasisj,onesfd,0,0,rangeval,wtfdj);
            Dmat(indexj) = Dmatj;
            %  loop through columns of CMAT
            mk2 = 0;
            for k=1:j
                betafdPark = betacell{k};
                if getestimate(betafdPark)
                    betafdk    = getfd(betafdPark);
                    betabasisk = getbasis(betafdk);
                    ncoefk     = length(getcoef(getfd(betafdPark)));
                    %  column indices of CMAT to fill
                    mk1 = mk2 + 1;
                    mk2 = mk2 + ncoefk;
                    indexk = mk1:mk2;
                    %  set up two weight functions
                    xfdk = xfdcell{k};
                    if wtconstant
                        xxfdjk = xfdj.*xfdk;
                    else
                        xxfdjk = (xfdj.*wt).*xfdk;
                    end
                    wtfdjk = sum(xxfdjk);
                    Cmatjk = inprod(betabasisj, betabasisk, 0, 0, rangeval, wtfdjk);
                    Cmat(indexj,indexk) = Cmatjk;
                    Cmat(indexk,indexj) = Cmatjk';
                end
            end
            %  attach penalty term to diagonal block
            lambda = getlambda(betafdParj);
            if lambda > 0
                Lfdj  = getLfd(betafdParj);
                Rmatj = eval_penalty(getbasis(getfd(betafdParj)), Lfdj);
                Cmat(indexj,indexj) = Cmat(indexj,indexj) + lambda.*Rmatj;
            end
        end
    end
    
    Cmat = (Cmat + Cmat')./2;
    
    %  check Cmat for singularity
    
    eigchk(Cmat);
    
    %  solve for coefficients defining BETA
    
    Cmatinv  = inv(Cmat);
    betacoef = Cmat\Dmat;
    
    %  set up fdPar objects for reg. fns. for BETAESTCELL
    
    betaestcell = betacell;
    mj2 = 0;
    for j=1:p
        betafdParj = betacell{j};
        if getestimate(betafdParj)
            betafdj = getfd(betafdParj);
            ncoefj  = size(getcoef(betafdj),1);
            mj1     = mj2 + 1;
            mj2     = mj2 + ncoefj;
            indexj  = mj1:mj2;
            betaestfdj     = putcoef(betafdj, betacoef(indexj));
            betaestfdPar   = putfd(betafdParj, betaestfdj);
        end
        betaestcell{j} = betaestfdPar;
    end
    
    %  set up fd object for predicted values in YHATFDOBJ
    
    nfine   = max(501,10*ynbasis+1);
    tfine   = linspace(rangeval(1), rangeval(2), nfine)';
    yhatmat = zeros(nfine,N);
    for j=1:p
        xmat    = eval_fd(tfine, xfdcell{j});
        betafdj = getfd(betaestcell{j});
        betavec = eval_fd(tfine, betafdj);
        yhatmat = yhatmat + xmat.*(betavec*ones(1,N));
    end
    yhatfd = smooth_basis(tfine, yhatmat, ybasisobj);
    
    df = NaN;
        
    if ~(isempty(y2cMap) || isempty(SigmaE))
    
        %  check dimensions of y2cMap and SigmaE

        y2cdim = size(y2cMap);
        if y2cdim(1) ~= ynbasis || ...
           y2cdim(2) ~= size(SigmaE, 1)  
             error('Dimensions of Y2CMAP not correct.');
        end
					
        ybasismat = eval_basis(tfine, ybasisobj);

        deltat    = tfine(2) - tfine(1);

        %  compute BASISPRODMAT

        basisprodmat = zeros(ncoef,ynbasis*N);

        mj2 = 0;
        for j = 1:p
            betafdParj = betacell{j};
            betabasisj = getbasis(getfd(betafdParj));
            ncoefj     = getnbasis(betabasisj);
            bbasismatj = eval_basis(tfine, betabasisj);
            xfdj       = xfdcell{j};
            tempj      = eval_fd(tfine, xfdj);
            %  row indices of BASISPRODMAT to fill
            mj1    = mj2 + 1;
            mj2    = mj2 + ncoefj;
            indexj = mj1:mj2;
            %  inner products of beta basis and response basis
            %    weighted by covariate basis functions
            mk2 = 0;
            for k = 1:ynbasis
                %  row indices of BASISPRODMAT to fill
                mk1    = mk2 + 1;
                mk2    = mk2 + N;
                indexk = mk1:mk2;
                tempk  = bbasismatj*ybasismat(:,k);
                basisprodmat(indexj,indexk) = ...
                                 deltat*crossprod(tempk,tempj);
            end
        end

        %  compute variances of regression coefficient function values
		
        c2bMap    = Cmat\basisprodmat;
        VarCoef   = y2cMap * SigmaE * y2cMap';
        CVariance = kron(VarCoef,diag(ones(N,1)));
        bvar      = c2bMap * CVariance * c2bMap';
        betastderrcell = cell(p,1);
        mj2 = 0;
        for j = 1:p
            betafdParj = betacell{j};
            betabasisj = getbasis(getfd(betafdParj));
            ncoefj     = getnbasis(betabasisj);
            mj1 	   = mj2 + 1;
            mj2 	   = mj2 + ncoefj;
            indexj 	   = mj1:mj2;
            bbasismat  = eval_basis(tfine, betabasisj);
            bvarj      = bvar(indexj,indexj);
            bstderrj   = sqrt(diag(bbasismat * bvarj * bbasismat'));
            bstderrfdj = smooth_basis(tfine, bstderrj, betabasisj);
            betastderrcell{j} = bstderrfdj;
        end
    else
        betastderrcell = [];
        bvar           = [];
        c2bMap         = [];
    end
    
else  
    
    %  ----------------------------------------------------------------
    %                   YFDPAR is scalar or multivariate
    %  ----------------------------------------------------------------
    
    ymat = yfdPar;
    
    Zmat  = [];
    Rmat  = [];
    pjsum = 0;
    pjvec = zeros(p,1);
    for j=1:p
        xfdj       = xfdcell{j};
        xcoef      = getcoef(xfdj);
        xbasis     = getbasis(xfdj);
        betafdParj = betacell{j};
        bbasis     = getbasis(getfd(betafdParj));
        bnbasis    = getnbasis(bbasis);
        pjvec(j)   = bnbasis;
        Jpsithetaj = inprod_basis(xbasis,bbasis);
        Zmat       = [Zmat,xcoef'*Jpsithetaj];
        if getestimate(betafdParj)
            lambdaj    = getlambda(betafdParj);
            if lambdaj > 0
                Lfdj  = getLfd(betafdParj);
                Rmatj = lambdaj.*eval_penalty(betafdParj, Lfdj);
            else
                Rmatj = zeros(pjvec(j));
            end
            Rmat  = [ [Rmat,    zeros(pjsum,bnbasis)]; ...
                      [zeros(bnbasis,pjsum), Rmatj] ];
            pjsum = pjsum + bnbasis;
        end
    end
    
    %  -----------------------------------------------------------
    %          set up the linear equations for the solution
    %  -----------------------------------------------------------
    
    %  solve for coefficients defining BETA
    
    if any(wt ~= 1)
        rtwt   = sqrt(wt);
        Zmatwt = Zmat.*(rtwt*ones(1,p));
        Cmat   = Zmatwt'*Zmatwt + Rmat;
        Dmat   = Zmatwt'*ymat;
    else
        Cmat = Zmat'*Zmat + Rmat;
        Dmat = Zmat'*ymat;
    end
    
    eigchk(Cmat);
    
    betacoef = Cmat\Dmat;
    
    % compute degrees of freedom measure
    
    Cmatinv = inv(Cmat);
    temp    = Cmat\Zmat';
    df      = sum(diag(Zmat*temp));
    
    %  set up fdPar object for BETAESTFDPAR
    
    betaestcell = betacell;
    mj2 = 0;
    for j=1:p
        mj1    = mj2 + 1;
        mj2    = mj2 + pjvec(j);
        indexj = mj1:mj2;
        betacoefj      = betacoef(indexj);
        betafdParj     = betacell{j};
        betafdj        = getfd(betafdParj);
        betaestfdj     = putcoef(betafdj, betacoefj);
        betaestfdParj  = putfd(betafdParj, betaestfdj);
        betaestcell{j} = betaestfdParj;
    end
    
    %  set up fd object for predicted values
    
    yhatmat = zeros(N,1);
    for j=1:p
        xfdj = xfdcell{j};
        if isa_fd(xfdj)    
            xbasis  = getbasis(xfdj);
            xnbasis = getnbasis(xbasis);
            xrng    = getbasisrange(xbasis);
            nfine   = max(501,10*xnbasis+1);
            tfine   = linspace(xrng(1), xrng(2), nfine)';
            deltat  = tfine(2)-tfine(1);
            xmat    = eval_fd(tfine, xfdj);
            betafdj = getfd(betaestcell{j});
            betamat = eval_fd(tfine, betafdj);
            yhatmat = yhatmat + deltat.*(xmat'*betamat    - ...
                      0.5.*(xmat(1,    :)'*betamat(1)     + ...
                            xmat(nfine,:)'*betamat(nfine)));
        else
            betaj   = getcoef(getfd(betaestcell{j}));
            yhatmat = yhatmat + xfdj*betaj';
        end
    end
    yhatfd = yhatmat;
        
    %  -----------------------------------------------------------------------
    %        Compute pointwise standard errors of regression coefficients
    %               if both y2cMap and SigmaE are supplied.
    %  -----------------------------------------------------------------------

    if ~(isempty(y2cMap) || isempty(SigmaE))
    
        %  compute linear mapping c2bMap takinging coefficients for
        %  response into coefficients for regression functions

        c2bMap = Cmat\Zmat';
        y2bmap = c2bMap;
        bvar   = y2bmap*SigmaE*y2bmap';
        betastderrcell = cell(p,1);
        mj2 = 0;
        for j = 1:p
	        betafdParj = betacell{j};
            betabasisj = getbasis(getfd(betafdParj));
	        ncoefj     = getnbasis(betabasisj);
            mj1        = mj2 + 1;
            mj2        = mj2 + ncoefj;
            indexj     = mj1:mj2;
            bvarj      = bvar(indexj,indexj);
            betarng    = getbasisrange(betabasisj);
            nfine      = max([501,10*ncoefj+1]);
            tfine      = linspace(betarng(1), betarng(2), nfine);
            bbasismat  = eval_basis(tfine, betabasisj);
            bstderrj   = sqrt(diag(bbasismat*bvarj*bbasismat'));
            bstderrfdj = smooth_basis(tfine, bstderrj, betabasisj);
            betastderrcell{j} = bstderrfdj;
        end
    else
        betastderrcell = [];
        bvar           = [];
        c2bMap         = [];
    end
end

%  -----------------------------------------------------------------------
%                  Set up output struct object
%  -----------------------------------------------------------------------

fRegressStruct.yfdPar     = yfdPar;
fRegressStruct.xfdcell    = xfdcell;
fRegressStruct.betacell   = betacell;
fRegressStruct.betahat    = betaestcell;
fRegressStruct.yhat       = yhatfd;
fRegressStruct.Cmat       = Cmat;
fRegressStruct.Dmat       = Dmat;
fRegressStruct.Cmatinv    = Cmatinv;
fRegressStruct.w          = wt;
fRegressStruct.df         = df;
fRegressStruct.y2cMap     = y2cMap;
fRegressStruct.SigmaE     = SigmaE;
fRegressStruct.betastderr = betastderrcell;
fRegressStruct.bvar       = bvar;
fRegressStruct.c2bMap     = c2bMap;

%  ------------------------------------------------------------

function eigchk(Cmat)
eigval = sort(eig(Cmat));
ncoef  = length(eigval);
if (eigval(1) < 0)
    disp('Smallest 10 eigenvalues:')
    neig = min([length(eigval),10]);
    for ieig=1:neig
        fprintf('  %g \n',eigval(ieig));
    end
    disp('Largest  10 eigenvalues:')
    for ieig=ncoef-neig+1:ncoef
        fprintf('  %g \n',eigval(ieig));
    end
    error('Negative eigenvalue of coefficient matrix.');
end
if (eigval(1) == 0)
    error('Zero eigenvalue of coefficient matrix.');
end
logcondition = log10(eigval(ncoef)) - log10(eigval(1));
if logcondition > 12
    warning('Wid1:singular', ...
        ['Near singularity in coefficient matrix.\n', ...
        'Eigenvalues range from ',num2str(eigval(1)), ...
        ' to ',num2str(eigval(ncoef)),'.']);
end

