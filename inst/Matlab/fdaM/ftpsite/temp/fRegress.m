function fRegressCell = ...
            fRegress(yfdPar, xfdcell, betacell)
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
%  
%  Returns:  
%  FREGRESSCELL   ...  A cell containing members with names:
%    yfdPar      ... first  argument of FREGRESS
%    xfdcell     ... second argument of FREGRESS
%    betacell    ... third  argument of FREGRESS
%    betaestcell ... estimated regression functions 
%    yhatfdobj   ... functional data object containing fitted functions 
%    Cmatinv     ... inverse of the coefficient matrix, needed for
%                    function FREGRESS.STDERR that computes standard errors

%  Last modified 20 July 2006

if nargin < 3
    error('Less than three arguments supplied.');
end

%  get number of independent variables 
    
p = length(xfdcell);

%  Check BETACELL
    
if length(betacell) ~= p
    error(['Number of regression coefficients does not match', ...
           ' number of independent variables.']);
end

for j=1:p
    betafdParj = betacell{j};
    if strcmp(class(betafdParj),'fd') 
        betafdParj  = fdPar(betafdParj);
        betacell{j} = betafdParj;
    end
    if ~strcmp(class(betafdParj),'fdPar')
        error(['BETACELL{',num2str(j),'} is not a FDPAR object.'])
    end
end

%  Get sample size and check YFDPAR.  

if strcmp(class(yfdPar), 'fdpar') || ...
   strcmp(class(yfdPar), 'fd')
    %  ----------------------------------------------------------------
    %                   YFDPAR is functional
    %  ----------------------------------------------------------------
    
    if strcmp(class(yfdPar), 'fd')
        yfdPar = fdPar(yfdPar);
    end
    yfdobj  = getfd(yfdPar);
    ycoef   = getcoef(yfdobj);
    if length(size(ycoef)) > 2
        error('YFDOBJ from YFDPAR is not univariate.');
    end
    N         = size(ycoef,2);
    ybasisobj = getbasis(yfdobj);
    rangeval  = getbasisrange(ybasisobj);
    ynbasis   = getnbasis(ybasisobj);
    nfine     = max(501,10*ynbasis+1);
    tfine     = linspace(rangeval(1), rangeval(2), nfine)';
    deltat    = tfine(2) - tfine(1);
    ymat      = eval_fd(tfine, yfdobj);
    
    %  get number of independent variables 
    
    p = length(xfdcell);
    
    %  check each cell.  If the object is a vector of length N,
    %  it is converted to a functional data object with a 
    %  constant basis
    
    onebasis = create_constant_basis(rangeval);
    
    for j=1:p
        xfdj = xfdcell{j};
        if isa_fd(xfdj)
            xcoef = getcoef(xfdj);
            if length(size(xcoef)) > 2
                error(['Covariate ',num2str(j),' is not univariate.']);
            end
            rangevalx  = getbasisrange(getbasis(xfdj));
            if any(rangevalx ~= rangeval)
                error(['Range for covariate ',num2str(j), ...
                        ' does not match that of YFDOBJ.']);
            end
        elseif strcmp(class(xfdj), 'double')
            xfdcell{j} = fd(xfdj(:)', onebasis);
        else
            error(['Covariate ', num2str(j),       ...
                   ' is neither a functional nor', ...
                   ' a multivariate object.']);
        end
        %  check size of coefficient array
        coefj = getcoef(xfdcell{j});
        Nj = size(coefj, 2);
        if Nj ~= N
            error('Incorrect number of replications in XFDCELL');
        end
    end
    
    if length(betacell) ~= p
        error(['Number of regression coefficients does not match', ...
                ' number of independent variables.']);
    end
    
    %  set up a matrix of values of covariates over a fine mesh
    
    xmat = zeros(nfine, N, p);
    betamatcell = cell(p,1);
    for j=1:p
        xmatj          = eval_fd(tfine, xfdcell{j});
        xmat(:,:,j)    = xmatj;
        betabasisj     = getbasis(getfd(betacell{j}));
        betamatcell{j} = eval_basis(tfine, betabasisj);
    end
    
    %  -----------------------------------------------------------
    %          set up the linear equations for the solution
    %  -----------------------------------------------------------
    
    %  compute the total number of coefficients to be estimated
    
    ncoef = 0;
    for j=1:p
        betafdParj = betacell{j};
        betafdj    = getfd(betafdParj);
        ncoefj     = size(getcoef(betafdj),1);
        ncoef      = ncoef + ncoefj;
    end
    
    Cmat = zeros(ncoef,ncoef);
    Dmat = zeros(ncoef,1);
    
    %  loop through rows of CMAT
    
    mj2 = 0;
    for j=1:p
        betafdParj = betacell{j};
        ncoefj     = length(getcoef(getfd(betafdParj)));
        %  row indices of CMAT and DMAT to fill
        mj1    = mj2 + 1;
        mj2    = mj2 + ncoefj;
        indexj = mj1:mj2;
        %  compute right side of equation DMAT
        xywtvec  = sum(squeeze(xmat(:,:,j)).*ymat,2);
        betamatj = betamatcell{j};
        temp     = betamatj.*(xywtvec*ones(1,ncoefj));
        Dmatj    = deltat.*sum(temp)';
        Dmat(indexj) = Dmatj;
        %  loop through columns of CMAT
        mk2 = 0;
        for k=1:j
            betafdPark = betacell{k};
            ncoefk     = length(getcoef(getfd(betafdPark)));
            %  column indices of CMAT to fill
            mk1 = mk2 + 1;
            mk2 = mk2 + ncoefk;
            indexk = mk1:mk2;
            %  set up two weight functions
            xxwtvec  = sum(xmat(:,:,j).*xmat(:,:,k),2);
            betamatk = betamatcell{k};
            temp     = betamatj.*(xxwtvec*ones(1,ncoefj));
            Cmatjk   = deltat.*temp'*betamatk;
            Cmat(indexj,indexk) = Cmatjk;
            Cmat(indexk,indexj) = Cmatjk';
        end
        %  attach penalty term to diagonal block
        lambda = getlambda(betafdParj);
        if lambda > 0
            Rmatj = eval_penalty(getbasis(getfd(betafdParj)));
            Cmat(indexj,indexj) = Cmat(indexj,indexj) + lambda.*Rmatj;
        end
    end
    
    Cmat = (Cmat + Cmat')./2;
    
    %  check Cmat for singularity
    
    eigval = sort(eig(Cmat));
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
        warning('Near singularity in coefficient matrix.');
        warning(['Eigenvalues range from ',num2str(eigval(1)), ...
                 ' to ',num2str(eigval(ncoef)),'.']);
    end
    
    %  solve for coefficients defining BETA
    
    Cmatinv  = inv(Cmat);
    betacoef = Cmatinv*Dmat;
    
    %  set up fdPar object for BETAFDPAR
    
    betaestcell = betacell;
    mj2 = 0;
    for j=1:p
        betafdParj     = betacell{j};
        betafdj        = getfd(betafdParj);
        ncoefj = size(getcoef(betafdj),1);
        mj1    = mj2 + 1;
        mj2    = mj2 + ncoefj;
        indexj = mj1:mj2;
        betaestfdj     = putcoef(betafdj, betacoef(indexj));
        betaestfdPar   = putfd(betafdParj, betaestfdj);
        betaestcell{j} = betaestfdPar;
    end
    
    %  set up fd object for predicted values
    
    yhatmat = zeros(nfine,N);
    for j=1:p
        xmat    = eval_fd(tfine, xfdcell{j});
        betafdj = getfd(betaestcell{j});
        betavec = eval_fd(tfine, betafdj);
        yhatmat = yhatmat + xmat.*(betavec*ones(1,N));
    end
    yhatfdobj = data2fd(yhatmat, tfine, ybasisobj);
        
elseif strcmp(class(yfdPar),'double')
    
    %  ----------------------------------------------------------------
    %                   YFDPAR is scalar or multivariate
    %  ----------------------------------------------------------------
    
    ymat = yfdPar;
    N    = size(ymat,1);
    
    %  check each cell.  If the object is a functional data object,
    %  it is converted to a multivariate object
    
    Zmat  = [];
    Rmat  = [];
    pjsum = 0;
    pjvec = zeros(p,1);
    for j=1:p
        xfdj = xfdcell{j};
        if strcmp(class(xfdj), 'fd')
            xcoef  = getcoef(xfdj);
            Nj     = size(xcoef,2);
            if Nj ~= N
                error(['Coefficient matrix ',num2str(j), ...
                       ' has the wrong number of columns.']);
            end
            xbasis     = getbasis(xfdj);
            betafdParj = betacell{j};
            bbasis     = getbasis(getfd(betafdParj));
            bnbasis    = getnbasis(bbasis);
            pjvec(j)   = bnbasis;
            Jpsithetaj = inprod_basis(xbasis,bbasis);
            Zmat       = [Zmat,xcoef'*Jpsithetaj];
            lambdaj    = getlambda(betafdParj);
            if getestimate(betafdParj) && lambdaj > 0
                Lfdj  = getLfd(betafdParj);
                Rmatj = lambdaj.*eval_penalty(betafdParj, Lfdj);
            else
                Rmatj = zeros(pjvec(j));
            end
            Rmat  = [ [Rmat,    zeros(pjsum,bnbasis)]; ...
                      [zeros(bnbasis,pjsum), Rmatj] ];
            pjsum = pjsum + bnbasis;
        elseif strcmp(class(xfdj), 'double')
            Zmatj    = xfdj;
            [Nj,pj]  = size(Zmatj);
            pjvec(j) = pj;
            if Nj ~= N
                error(['Covariate matrix ',num2str(j), ...
                       ' has the wrong number of rows.']);
            end
            Zmat  = [Zmat,Zmatj];
            Rmatj = zeros(pj);
            Rmat  = [ [Rmat,  zeros(pjsum,pj)]; ...
                      [zeros(pj,pjsum), Rmatj] ];
            pjsum = pjsum + pj;
        else
            error(['Covariate ',num2str(j), ...
                   ' is neither a functional nor a multivariate object.']);
        end
    end
    
    %  -----------------------------------------------------------
    %          set up the linear equations for the solution
    %  -----------------------------------------------------------
    
    %  solve for coefficients defining BETA
    
    Cmat     = Zmat'*Zmat + Rmat;
    Cmatinv  = inv(Cmat);
    Dmat     = Zmat'*ymat;
    
    betacoef = Cmatinv*Dmat;
    
    %  set up fdPar object for BETAESTFDPAR
    
    betaestcell = betacell;
    onebasis    = create_constant_basis([0,1]);
    mj2 = 0;
    for j=1:p
        mj1 = mj2 + 1;
        mj2 = mj2 + pjvec(j);
        indexj = mj1:mj2;
        betacoefj  = betacoef(indexj);
        betafdParj = betacell{j};
        if strcmp(class(xfdj), 'fd')        
            betafdj        = getfd(betafdParj);
            betaestfdj     = putcoef(betafdj, betacoefj);
            betaestfdParj  = putfd(betafdParj, betaestfdj);
            betaestcell{j} = betaestfdParj;
        else
            betaestfdj     = fd(betacoefj',onebasis);
            betaestfdParj  = putfd(betafdParj, betaestfdj);
            betaestcell{j} = betaestfdParj;
        end
    end
    
    %  set up fd object for predicted values
    
    yhatmat = zeros(N,1);
    for j=1:p
        xfdj = xfdcell{j};
        if strcmp(class(xfdj), 'fd')    
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
    yhatfdobj = yhatmat;
        
else
    %  YFDOBJ is neither functional nor multivariate
    error('YFDOBJ is neither functional nor multivariate.');
end

fRegressCell{1} = yfdPar;
fRegressCell{2} = xfdcell;
fRegressCell{3} = betacell;
fRegressCell{4} = betaestcell;
fRegressCell{5} = yhatfdobj;
fRegressCell{6} = Cmatinv;
