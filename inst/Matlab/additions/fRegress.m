function fRegressCell = ...
            fRegress(yfdPar, xfdcell, betacell, wt)
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
%  
%  Returns:  
%  FREGRESSCELL  ...  A cell containing members with names:
%    yfdPar      ... first  argument of FREGRESS
%    xfdcell     ... second argument of FREGRESS
%    betacell    ... third  argument of FREGRESS
%    betaestcell ... estimated regression functions 
%    yhatfdobj   ... functional data object containing fitted functions 
%    Cmatinv     ... inverse of the coefficient matrix, needed for
%                    function FREGRESS_STDERR that computes standard errors
%    wt          ... weights for observations
%    df          ... degrees of freedom for fit  (scalar response)
%    OCV         ... ordinary cross validation score (scalar response)
%    GCV         ... generalized cross validation score (scalar resonse)

%  Last modified 15 May 2009

if nargin < 3
    error('Less than three arguments supplied.');
end

%  Check that xfdcell is a cell object

if ~iscell(xfdcell)
    error('Argument XFDCELL is not a cell object.');
end

%  get number of independent variables 
    
p = length(xfdcell);

%  check contents of XFDCELL

xerror = 0;
for j=1:p
    xfdj = xfdcell{j};
    if ~(isa_fd(xfdj) || isa_fdPar(xfdj) || isnumeric(xfdj))
        disp(['XFDCELL{',num2str(j), ...
               '} is not a vector, an FD object or an FDPAR object.']);
        xerror = 1;
    end
end

%  Check BETACELL
  
if ~iscell(betacell)
    error('Argument BETACELL is not a cell object.');
end

if length(betacell) ~= p
    error(['Number of regression coefficients does not match', ...
           ' number of independent variables.']);
end

%  check contents of BETAFDCELL

berror = 0;
for j=1:p
    betafdParj = betacell{j};
    if isa_fd(betafdParj) 
        betafdParj  = fdPar(betafdParj);
        betacell{j} = betafdParj;
    end
    if ~isa_fdPar(betafdParj)
        disp(['BETACELL{',num2str(j),'} is not a FDPAR object.']);
        berror = 1;
    end
end

if xerror || berror
    error('An error has been found in either XFDCELL or BETACELL.'); 
end

%  set up a constant basis for vector independent variables

betafdPar = betacell{1};
betafd    = getfd(betafdPar);
betabasis = getbasis(betafd);
betarange = getbasisrange(betabasis);
onebasis  = create_constant_basis(betarange);
onesfd    = fd(1,onebasis);

%  Get sample size and check YFDPAR.  

if isa_fdPar(yfdPar) || isa_fd(yfdPar)
    
    %  ----------------------------------------------------------------
    %                   YFDPAR is functional
    %  ----------------------------------------------------------------
    
    if isa_fd(yfdPar)
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
    
    %  check weight vector WT
    
    if nargin < 4
        wt = ones(N,1);
    end
    
    if length(wt) ~= N
        error('The number of observation weights is incorrect.');
    end
    
    if any(wt < 0)
        error('Negative observation weights encountered.');
    end
    
    if var(wt) > 0 
        wtconstant = 0;
    else
        wtconstant = 1;
    end
        
    %  check each cell.  If the object is a vector of length N,
    %  it is converted to a functional data object with a 
    %  constant basis
    
    xerror = 0;
    for j=1:p
        xfdj = xfdcell{j};
        if isa_fd(xfdj)
            xcoef = getcoef(xfdj);
            if length(size(xcoef)) > 2
                error(['Covariate ',num2str(j),' is not univariate.']);
            end
            rangevalx  = getbasisrange(getbasis(xfdj));
            if any(rangevalx ~= rangeval)
                disp(['Range for covariate ',num2str(j), ...
                        ' does not match that of YFDOBJ.']);
                xerror = 1;
            end
        elseif strcmp(class(xfdj), 'double')
            xfdcell{j} = fd(xfdj(:)', onebasis);
        else
            disp(['Covariate ', num2str(j),       ...
                   ' is neither a functional nor', ...
                   ' a multivariate object.']);
            xerror = 1;
        end
        %  check size of coefficient array
        coefj = getcoef(xfdcell{j});
        Nj = size(coefj, 2);
        if Nj ~= N
            disp('Incorrect number of replications in XFDCELL');
            xerror = 1;
        end
    end
    if xerror, error(''); end
    
    if length(betacell) ~= p
        error(['Number of regression coefficients does not match', ...
                ' number of independent variables.']);
    end
    
    %  check weights

    if length(wt) ~= N 
        error('Number of weights not equal to N.');
    end
    if any(wt < 0)    
        error('Negative weights found.');
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
                Rmatj = eval_penalty(getbasis(getfd(betafdParj)));
                Cmat(indexj,indexj) = Cmat(indexj,indexj) + lambda.*Rmatj;
            end
        end
    end
    
    Cmat = (Cmat + Cmat')./2;
    
    %  check Cmat for singularity
    
    eigchk(Cmat);
    
    %  solve for coefficients defining BETA
    
    Cmatinv  = inv(Cmat);
    betacoef = Cmatinv*Dmat;
    
    %  set up fdPar object for BETAFDPAR
    
    betaestcell = betacell;
    mj2 = 0;
    for j=1:p
        betafdParj     = betacell{j};
        if getestimate(betafdParj)
            betafdj        = getfd(betafdParj);
            ncoefj = size(getcoef(betafdj),1);
            mj1    = mj2 + 1;
            mj2    = mj2 + ncoefj;
            indexj = mj1:mj2;
            betaestfdj     = putcoef(betafdj, betacoef(indexj));
            betaestfdPar   = putfd(betafdParj, betaestfdj);
        end
        betaestcell{j} = betaestfdPar;
    end
    
    %  set up fd object for predicted values
    
    nfine     = max(501,10*ynbasis+1);
    tfine     = linspace(rangeval(1), rangeval(2), nfine)';
    yhatmat = zeros(nfine,N);
    for j=1:p
        xmat    = eval_fd(tfine, xfdcell{j});
        betafdj = getfd(betaestcell{j});
        betavec = eval_fd(tfine, betafdj);
        yhatmat = yhatmat + xmat.*(betavec*ones(1,N));
    end
    yhatfdobj = data2fd(yhatmat, tfine, ybasisobj);
    
    df = NaN;
    OCV = NaN;
    GCV = NaN;
        
elseif strcmp(class(yfdPar),'double')
    
    %  ----------------------------------------------------------------
    %                   YFDPAR is scalar or multivariate
    %  ----------------------------------------------------------------
    
    ymat = yfdPar;
    N    = size(ymat,1);
    
    %  check weight vector WT
    
    if nargin < 4
        wt = ones(N,1);
    end
    
    if length(wt) ~= N
        error('The number of observation weights is incorrect.');
    end
    
    if any(wt < 0)
        error('Negative observation weights encountered.');
    end

    %  check each cell.  If the object is a functional data object,
    %  it is converted to a multivariate object
    
    Zmat  = [];
    Rmat  = [];
    pjsum = 0;
    pjvec = zeros(p,1);
    for j=1:p
        xfdj = xfdcell{j};
        if isa_fd(xfdj)
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
            if getestimate(betafdParj)
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
            end
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
    
    if any(wt ~= 1)
        rtwt   = sqrt(wt);
        Zmatwt = Zmat.*(rtwt*ones(1,p));
        Cmat   = Zmatwt'*Zmatwt + Rmat;
        Dmat   = Zmatwt'*ymatwt;
    else
        Cmat = Zmat'*Zmat + Rmat;
        Dmat = Zmat'*ymat;
    end
    
    eigchk(Cmat);
    
    Cmatinv  = inv(Cmat);
    betacoef = Cmatinv*Dmat;
    
    % compute degrees of freedom measure
    
%    df = sum(diag(Zmat*Cmatinv*Zmat'));
    hatvals = diag(Zmat * Cmatinv * Zmat');
    df = sum(hatvals);
    
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
        if isa_fd(xfdj)        
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
    yhatfdobj = yhatmat;

    OCV = sum( (ymat-yhatmat).^2./(1-hatvals).^2 );
    GCV = sum( (ymat-yhatmat).^2 )/( (sum(1-hatvals)).^2 );
    
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
fRegressCell{7} = wt;
fRegressCell{8} = df;
fRegressCell{9} = OCV;
fRegressCell{10} = GCV;

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

