function stderrStruct = fRegress_stderr(fRegressStruct, y2cMap, SigmaE)

%  FREGRESS.STDERR  computes standard error estimates for regression
%       coefficient functions estimated by function FREGRESS.
%
%  Arguments:
%  
%  FREGRESSCELL ... a cell object produced by function FREGRESS.
%  Y2CMAP       ... the matrix mapping from the vector of observed values
%                   to the coefficients for the dependent variable.  
%                   This is output by function SMOOTH_BASIS.  If this is 
%                   supplied, confidence limits are computed, otherwise not.
%  SIGMAE       ... Estimate of the covariances among the residuals.  This
%                   can only be estimated after a preliminary analysis 
%                   with FREGRESS.
%  
%  Returns:  
%  BETASTDERRCELL ... a cell array object, each element containing a fdPar 
%                     object for the standard error of a regression function.
%  BVARIANCE      ... the symmetric matrix of sampling variances and
%                     covariances for the matrix of regression coefficients
%                     for the regression functions.  These are stored
%                     column-wise in defining BVARIANCE.  This is only 
%                     computed if the argument Y2CMAP is supplied.
%  C2BMAP         ... the matrix mapping from response variable coefficients 
%                     to coefficients for regression coefficients

%  Last modified 20 July 2006

%  get number of independent variables 
    
yfdPar   = fRegressStruct.yfdPar;
xfdcell  = fRegressStruct.xfdcell;
betacell = fRegressStruct.betacell;
Cmatinv  = fRegressStruct.Cmatinv;

p = length(xfdcell);

%  compute number of coefficients

ncoef = 0;
for j = 1:p
    betaParfdj = betacell{j};
    ncoefj     = getnbasis(getbasis(getfd(betaParfdj)));
    ncoef      = ncoef + ncoefj;
end

if strcmp(class(yfdPar), 'fdPar') || strcmp(class(yfdPar), 'fd')

    %  ----------------------------------------------------------------
    %           YFDPAR is functional for a functional parameter
    %  ----------------------------------------------------------------

    if strcmp(class(yfdPar), 'fd'), yfdPar = fdPar(yfdPar);  end

    %  get number of replications and basis information for YFDPAR

    yfd       = getfd(yfdPar);
    N         = size(getcoef(yfd),2);
    ybasisobj = getbasis(getfd(yfdPar));
    rangeval  = getbasisrange(ybasisobj);
    ynbasis   = getnbasis(ybasisobj);
    nfine     = max([501,10*ynbasis+1]);
    tfine     = linspace(rangeval(1), rangeval(2), nfine)';
    ybasismat = eval_basis(tfine, ybasisobj);
    deltat    = tfine(2) - tfine(1);

    %  compute BASISPRODMAT

    basisprodmat = zeros(ncoef,ynbasis*N);

    mj2 = 0;
    for j=1:p
        betafdParj = betacell{j};
        betafdj    = getfd(betafdParj);
        betabasisj = getbasis(betafdj);
        ncoefj     = getnbasis(betabasisj);
        bbasismat  = eval_basis(tfine, betabasisj);
        xfdj       = xfdcell{j};
        xmatj      = eval_fd(tfine,xfdj);
        %  row indices of BASISPRODMAT to fill
        mj1    = mj2 + 1;
        mj2    = mj2 + ncoefj;
        indexj = mj1:mj2;
        %  inner products of beta basis and response basis
        %    weighted by covariate basis functions
        mk2 = 0;
        for k=1:ynbasis
            mk1    = mk2 + 1;
            mk2    = mk2 + N;
            indexk = mk1:mk2;
            wtmatk = ybasismat(:,k)*ones(1,ncoefj);
            tempk  = bbasismat .* wtmatk;
            basisprodmat(indexj,indexk) = deltat.*tempk'*xmatj;
        end
    end

    %  check dimensions of Y2CMAP

    y2cdim = size(y2cMap);
    if y2cdim(1) ~= ynbasis || y2cdim(2) ~= size(SigmaE,1)
        error('Dimensions of Y2CMAP not correct.');
    end

    %  compute variances of regression coefficient function values

    c2bMap    = Cmatinv*basisprodmat;
    VarCoef   = y2cMap*SigmaE*y2cMap';
    CVariance = kron(VarCoef,eye(N));
    bvar = c2bMap*CVariance*c2bMap';
    mj2 = 0;
    betastderrcell = cell(p,1);
    for j=1:p
        betafdParj = betacell{j};
        betafdj    = getfd(betafdParj);
        betabasisj = getbasis(betafdj);
        ncoefj     = getnbasis(betabasisj);
        mj1 = mj2 + 1;
        mj2 = mj2 + ncoefj;
        indexj = mj1:mj2;
        bbasismat  = eval_basis(tfine, betabasisj);
        bvarj      = bvar(indexj,indexj);
        bstderrj   = sqrt(diag(bbasismat*bvarj*bbasismat'));
        bstderrfdj = data2fd(bstderrj, tfine, betabasisj);
        betastderrcell{j} = bstderrfdj;
    end

else
    %  ----------------------------------------------------------------
    %                   YFDPAR is scalar or multivariate
    %  ----------------------------------------------------------------

    Zmat = [];
    for j = 1:p
        xfdj = xfdcell{j};
        if strcmp(class(xfdj), 'fd')
            xcoef      = getcoef(xfdj);
            xbasis     = getbasis(xfdj);
            betafdParj = betacell{j};
            bbasis     = getbasis(getfd(betafdParj));
            Jpsithetaj = inprod(xbasis,bbasis);
            Zmat       = [Zmat, xcoef'*Jpsithetaj];
        else
            if strcmp(class(xfdj), 'double')
                Zmatj = xfdj;
                Zmat  = [Zmat,Zmatj];
            end
        end
    end

    %  compute linear mapping c2bMap takinging coefficients for
    %  response into coefficients for regression functions

    c2bMap = Cmatinv*Zmat';
    y2bmap = c2bMap;
    bvar   = y2bmap*SigmaE*y2bmap';
    mj2 = 0;
    for j=1:p
        betafdParj = betacell{j};
        betafdj    = getfd(betafdParj);
        ncoefj     = size(getcoef(betafdj),1);
        mj1 = mj2 + 1;
        mj2 = mj2 + ncoefj;
        indexj = mj1:mj2;
        bvarj = bvar(indexj,indexj);
        xfdj  = xfdcell{j};
        if strcmp(class(xfdj), 'fd')
            betabasisj = getbasis(betafdj);
            bnbasis    = getnbasis(betabasisj);
            betarng    = getbasisrange(betabasisj);
            nfine      = max(501,10*bnbasis+1);
            tfine      = linspace(betarng(1), betarng(2), nfine)';
            bbasismat  = eval_basis(tfine, betabasisj);
            bstderrj   = sqrt(diag(bbasismat*bvarj*bbasismat'));
            bstderrfdj = data2fd(bstderrj, tfine, betabasisj);
        else
            bstderrj   = sqrt(diag(bvarj));
            bstderrfdj = fd(bstderrj', onebasis);
        end
        betastderrcell{j} = bstderrfdj;
    end
end

stderrStruct.betastderr = betastderrcell;
stderrStruct.bvar       = bvar;
stderrStruct.c2bMap     = c2bMap;
