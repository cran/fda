function [Wfdobj, betm, yhatfd, Fstr, y2cMap, argvals, y] = ...
          smooth_monotone(argvals, y, fdParobj, zmat, wtvec, ...
                          conv, iterlim, active, dbglev)
%SMOOTH_MONOTONE smooths the relationship of Y to ARGVALS by fitting 
%     a monotone function of the form
%                   f(x) = b_0 + b_1 D^{-1} exp W(x)
%     where  W  is a function defined over the same range as ARGVALS, and
%             W + ln b_1 = log Df and w = D W = D^2f/Df.
%  The constant term b_0 in turn can be a linear combination of covariates:
%                       b_0 = zmat * c,
%  in which case coefficient vector C is returned in its place.
%  By default (ZMAT = []), it is just a constant.
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum w_i[y_i - f(x_i)]^2 +
%                     \lambda * \int [L W(x)]^2 dx
%  The function W(x) is expanded by the basis in functional data object
%    WFD.   The coefficients of this expansion are called "coefficients"
%    in the comments, while the b's are called "regression coefficients"
%
%  Arguments are ...
%  ARGVALS ...  Argument value array of length N, where N is the number of 
%               observed curve values for each curve.  It is assumed that
%               that these argument values are common to all observed
%               curves.  If this is not the case, you will need to 
%               run this function inside one or more loops, smoothing
%               each curve separately.
%  Y       ...  Function value array (the values to be fit).
%               If the functional data are univariate, this array will be 
%               an N by NCURVE matrix, where N is the number of observed
%               curve values for each curve and NCURVE is the number of
%               curves observed.
%               If the functional data are muliivariate, this array will be 
%               an N by NCURVE by NVAR matrix, where NVAR the number of
%               functions observed per case.  For example, for the gait
%               data, NVAR = 2, since we observe knee and hip angles.
%  FDPAROBJ ... A functional parameter or fdPar object.  This object 
%               contains the specifications for the functional data
%               object to be estimated by smoothing the data.  See
%               comment lines in function fdPar for details.
%               The functional data object WFD in FDPAROBJ is used
%               to initialize the optimization process.
%               Its coefficient array contains the starting values for 
%               the iterative minimization of mean squared error.
%  ZMAT    ...  a matrix of covariate values for the constant term.
%               It defaults to N one's if ZMAT is empty (ZMAT = []).
%  WTVEC   ...  a vector of weights, a vector of N one's by default.
%               WTVEC may also be an order N matrix (usually symmetric)
%  CONV    ...  convergence criterion, 0.0001 by default
%  ITERLIM ...  maximum number of iterations, 50 by default.
%  ACTIVE  ...  indices among 1:NBASIS of parameters to optimize.
%               Defaults to 1:NBASIS.
%  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
%               no output, if 1, output at each iteration, if higher, 
%               output at each line search iteration. 1 by default.
%
%  Returns are:
%  WFD     ...  Functional data object for W. 
%               Its coefficient matrix an N by NCURVE (by NVAR) matrix
%               (or array), depending on whether the functional
%               observations are univariate or multivariate.
%  BETA    ...  The regression coefficients b_0 and b_1 for each
%               smoothed curve.  
%               If the curves are univariate and 
%                  ... ZMAT is empty, BETA is a 2 by NCURVE matrix.
%                  ... ZMAT has P columns, BETA is a P+2 by NCURVE matrix.
%               If the curves are multiivariate and 
%                  ... ZMAT is empty, BETA is a 2 by NCURVE by NVAR array.
%                  ... ZMAT has P columns, BETA is a P+2 by NCURVE by NVAR
%                      array.    
%  YHATFD  ...  A functional data object for the fitting function.
%               This is constructed using the basis for WFD, and this
%               basis may well be too simple to accommodate the curvature
%               in the monotone function that WFD defines.  It may be 
%               necessary to discard this object and use a richer basis
%               externally to smooth the values defined by
%                  beta(1) + beta(2).*eval_mon(evalarg, Wfd).
%  FSTR    ...  A struct object or a cell array of struct objects, one for 
%               each curve (and each variable if functions multivariate).
%               Each struct object has slots:
%                 f    ... The sum of squared errors   
%                 grad ... The gradient  
%                 norm ... The norm of the gradient  
%  Y2CMAP   ... For each estimated curve (and variable if functions are
%               multivariate, this is an N by NBASIS matrix containing
%               a linear mappping from data to coefficients used
%               for computing point-wise confidence intervals.  
%               If NCURVE = NVAR = 1, a matrix is returned.  Otherwise
%               an NCURVE by NVAR cell array is returned, with each 
%               entry being this mapping.
%  ARGVALS  ... Input argument values
%  Y        ... Input data to be smoothed

%  Last modified 1 September 2011
  
if nargin < 3
    error('The first three arguments are not specified.');
end

%  check ARGVALS

[argvals, n] = argcheck(argvals);
oneobs = ones(n,1); 

%  at least three points are necessary for monotone smoothing

if n < 3
    error('ARGVALS does not contain at least three values.');
end

%  check Y

[y, ncurve, nvar, ndim] = ycheck(y, n);

%  check fdParobj and get LAMBDA

fdParobj = fdParcheck(fdParobj);
lambda   = getlambda(fdParobj);

%  set up LFDOBJ

Lfdobj = getLfd(fdParobj);
Lfdobj = int2Lfd(Lfdobj);

%  set up BASIS

Wfd0     = getfd(fdParobj);
basisobj = getbasis(Wfd0);
nbasis   = getnbasis(basisobj);

%  set up initial coefficient array

coef0 = getcoef(Wfd0);  

%  initialize arguments that are not included

if nargin < 9, dbglev  = 1;        end
if nargin < 8, active  = 1:nbasis; end
if nargin < 7, iterlim = 50;       end
if nargin < 6, conv    = 0.0001;   end
if nargin < 5, wtvec   = oneobs;   end
if nargin < 4, zmat    = [];       end

%  check WTVEC

[wtvec, onewrd, matwrd] = wtcheck(n, wtvec);

%  check ZMAT and set NCOV

if ~isempty(zmat)
    zdim = size(zmat);
    if zdim(1) ~= n
        error('First dimension of ZMAT not correct.')
    end
    ncov = zdim(2); %  number of covariates
else
    ncov = 1;
end

%  set up some variables

ncovp1 = ncov + 1;
if ~matwrd
    wtroot = sqrt(wtvec);
end
climit = 100.0.*([-1; 1]*ones(1,nbasis));

%  set an indicator array for the fixed coefficients

inactive  = ones(1,nbasis);
inactive(active) = 0;
inactive  = find(inactive);

%  set up cell for storing basis function values

basiscell = cell(1,15);

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda*eval_penalty(basisobj, Lfdobj);
else
    Kmat  = zeros(nbasis,nbasis);
end

%  --------------------------------------------------------------------
%              loop through variables and curves
%  --------------------------------------------------------------------

%  set up arrays and cell arrays to contain returned information

if ndim == 2
    coef = zeros(nbasis,ncurve);
    betm = zeros(ncovp1,ncurve);
else
    coef = zeros(nbasis,ncurve,nvar);
    betm = zeros(ncovp1,ncurve,nvar);
end

if nargout > 2
    if ncurve > 1 || nvar > 1
        Fstr = cell(ncurve,nvar);
    else
        Fstr = [];
    end
end

if nargout > 4
    if ncurve > 1 || nvar > 1
        y2cMap = cell(ncurve,nvar);
    else
        y2cMap = [];
    end
end

for ivar=1:nvar
    for icurve=1:ncurve
        if ndim == 2
            yi    = squeeze(y(:,icurve));
            cveci = squeeze(coef0(:,icurve));
        else
            yi    = squeeze(y(:,icurve,ivar));
            cveci = squeeze(coef0(:,icurve,ivar));
        end

        %  Compute initial function and gradient values

        [Fstri, betmi, Dyhat, basiscell] = ...
            fngrad(yi, argvals, zmat, wtvec, cveci, lambda, ...
                   basisobj, Kmat, basiscell, inactive);

        %  compute the initial expected Hessian

        hessmat = hesscal(betmi, Dyhat, wtvec, lambda, Kmat, inactive);

        %  evaluate the initial line search direction vector

        deltac = linesearch(Fstri, hessmat, dbglev);

        %  initialize iteration status arrays

        iternum = 0;
        status = [iternum, Fstri.f, Fstri.norm, betmi'];
        if dbglev >= 1
            fprintf(['\nResults for curve ',num2str(icurve), ...
                   ' and variable ',num2str(ivar),'\n'])
            fprintf('\nIter.   PENSSE   Grad Length Intercept   Slope\n')
            fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
                [status(1:4),betmi(ncovp1)]);
        end
        if dbglev > 2
            for ibasis = 1:nbasis, fprintf('%10.4f%', cveci(ibasis)); end
            fprintf('\n');
            for ibetm  = 1:3,      fprintf('%10.4f%', betmi(ibetm)); end
            fprintf('\n');
        end

        %  ---------------------  Begin main iterations  ---------------

        MAXSTEPITER = 10;
        MAXSTEP = 100;
        trial   = 1;
        reset   = 0;
        linemat = zeros(3,5);
        betmold = betmi;
        cvecold = cveci;
        Foldstr = Fstri;
        dbgwrd  = dbglev >= 2;

        %  ---------------  beginning of optimization loop  -----------

        for iter = 1:iterlim
            iternum = iternum + 1;
            %  initialize logical variables controlling line search
            dblwrd = [0,0];  limwrd = [0,0];  ind = 0; ips = 0;
            %  compute slope at 0 for line search
            Fstri = Foldstr;
            linemat(2,1) = sum(deltac.*Fstri.grad);
            %  normalize search direction vector
            sdg          = sqrt(sum(deltac.^2));
            deltac       = deltac./sdg;
            linemat(2,1) = linemat(2,1)/sdg;
            % initialize line search vectors
            linemat(:,1:4) = [0; linemat(2,1); Fstri.f]*ones(1,4);
            stepiter  = 0;
            if dbglev >= 2
                fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                    [stepiter, linemat(:,1)']);
            end
            %  return with error condition if initial slope is nonnegative
            if linemat(2,1) >= 0
                if dbglev >= 2, disp('Initial slope nonnegative.'); end
                break;
            end
            %  return successfully if initial slope is very small
            if linemat(2,1) >= -1e-7;
                if dbglev >= 2, disp('Initial slope too small'); end
                break;
            end
            %  first step set to trial
            linemat(1,5)  = trial;
            %  -------  Begin line search iterations  -----------
            cvecnew = cveci;
            for stepiter = 1:MAXSTEPITER
                %  check that step size does not go beyond limits
                [linemat(1,5), ind, limwrd] = ...
                    stepchk(linemat(1,5), cveci, deltac, limwrd, ind, ...
                    climit, dbgwrd);
                if ind == 1, break; end  % limit break hit twice in a row
                if linemat(1,5) <= 1e-7
                    %  Current step size too small ... terminate
                    if dbglev >= 2
                        fprintf( ...
                            'Stepsize too small: %15.7f\n', linemat(1,5));
                    end
                end
                %  compute new function value and gradient
                cvecnew = cveci + linemat(1,5).*deltac;  
                [Fstri, betmi, Dyhat, basiscell] = ...
                    fngrad(yi, argvals, zmat, wtvec, cvecnew, lambda, ...
                           basisobj, Kmat, basiscell, inactive);
                linemat(3,5) = Fstri.f;
                %  compute new directional derivative
                linemat(2,5) = sum(deltac.*Fstri.grad);
                if dbglev >= 2
                    fprintf( ...
                        '                %3.f %10.4f %12.6f %12.6f\n', ...
                        [stepiter, linemat(:,5)']);
                end
                %  compute next line search step, test for convergence
                [linemat, ips, ind, dblwrd] = ...
                    stepit(linemat, ips, dblwrd, MAXSTEP);
                trial  = linemat(1,5);
                if trial == MAXSTEP, break; end
                %  ind == 0 means convergence
                if ind == 0 || ind == 5, break; end
            end
            %  -------  End of line search iterations  -----------
            cveci = cvecnew;
            %  check that function value has not increased
            if Fstri.f > Foldstr.f
                %  if it has, terminate iterations with a warning
                if dbglev >= 2
                    fprintf('Criterion increased:');
                    fprintf('%10.4f %10.4f\n',[Foldstr.f, Fstri.f]);
                end
                %  reset parameters and fit
                betmi  = betmold;
                cveci  = cvecold;
                Fstri  = Foldstr;
                deltac = -Fstri.grad;
                if dbglev > 2
                    for i = 1:nbasis, fprintf('%10.4f%', cveci(i)); end
                    fprintf('\n');
                    for i = 1:3,      fprintf('%10.4f%', betmi(i)); end
                    fprintf('\n');
                end
                if reset == 1
                    %  This is the second time in a row that this
                    %     has happened ...  quit
                    if dbglev >= 2
                        fprintf('Reset twice, terminating.\n');
                        return;
                    end
                else
                    reset = 1;
                end
            else
                %  function value has not increased,  check for convergence
                if abs(Foldstr.f-Fstri.f) < conv
                    status = [iternum, Fstri.f, Fstri.norm, betmi'];
                    if dbglev >= 1
                        fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
                            [status(1:4),betmi(ncovp1)]);
                    end
                    break;
                end
                %  update old parameter vectors and fit structure
                cvecold = cveci;
                betmold = betmi;
                Foldstr = Fstri;
                %  update the expected Hessian
                hessmat = ...
                    hesscal(betmi, Dyhat, wtvec, lambda, Kmat, inactive);
                %  update the line search direction vector
                deltac = linesearch(Fstri, hessmat, dbglev);
                reset = 0;
            end
            %  store iteration status
            status = [iternum, Fstri.f, Fstri.norm, betmi'];
            if dbglev >= 1
                fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
                    [status(1:4),betmi(ncovp1)]);
            end
        end
        
        %  save coefficients in arrays COEF and BETA
        
        if ndim == 2
            coef(:,icurve) = cveci;
            betm(:,icurve) = betmi;
        else
            coef(:,icurve,ivar) = cveci;
            betm(:,icurve,ivar) = betmi;
        end
        
        %  save Fstr if required in cell array.
        
        if nargout > 2
            if ncurve == 1 && nvar == 1
                Fstr = Fstri;
            else
                Fstr{icurve,ivar} = Fstri;
            end
        end
        
        %  save y2cMap if required in cell array.
        
        if nargout > 4
            y2cMapij = ((Dyhat' * Dyhat + lambda*Kmat)\Dyhat')./sqrt(n);
            if ncurve == 1 && nvar == 1
                y2cMap = y2cMapij;
            else
                y2cMap{icurve,ivar} = y2cMapij;
            end
        end
    end
end

Wfdobj = fd(coef, basisobj);

%  Set up yhatfd, a functional data object for the monotone curves
%  fitting the data.  
%  This can only be done if the covariate matrix ZMAT is NULL, meaning that
%  the same constant term is used for all curve values.

if isempty(zmat) 

  rangeval = getbasisrange(basisobj);
  narg     = 10*nbasis+1;
  evalarg  = linspace(rangeval(1), rangeval(2), narg)';
  hmat     = eval_mon(evalarg, Wfdobj);

  if ndim == 2 
    yhatmat = zeros(narg,ncurve);
    for icurve = 1:ncurve 
      yhatmat(:,icurve) = betm(1,icurve) + ...
                          betm(2,icurve)*hmat(:,icurve);
    end
    yhatcoef = project_basis(yhatmat, evalarg, basisobj);
    yhatfd   = fd(yhatcoef, basisobj);
  else 
    yhatcoef = zeros(nbasis,ncurve,nvar);
    yhatmati = zeros(narg,ncurve);
    for ivar = 1:nvar 
      for icurve = 1:ncurve 
        yhatmati(:,icurve) = betm(1,icurve,ivar) + ...
                             betm(2,icurve,ivar)*hmat(:,icurve,ivar);
      end
      yhatcoef(:,:,ivar) = project_basis(yhatmati, evalarg, basisobj);
    end
    yhatfd = fd(yhatcoef, basisobj);
  end 
else 
  yhatfd = [];
end

%  ----------------------------------------------------------------

    function [deltac, cosangle] = linesearch(Fstri, hessmat, dbglev)
        deltac   = -hessmat\Fstri.grad;
        cosangle = -sum(Fstri.grad.*deltac)./ ...
            sqrt(sum(Fstri.grad.^2).*sum(deltac.^2));
        if dbglev >= 2
            fprintf('Cos(angle) = %8.4f\n', cosangle);
        end
        if sum(cosangle) < 1e-7
            if dbglev >=2, fprintf('\nCosine of angle too small\n'); end
            deltac = -Fstri.grad;
        end
        
        %  ----------------------------------------------------------------
        
        function [Fstri, betmi, Dyhat, basiscell] = ...
                fngrad(yi, argvals, zmat, wtvec, cveci, lambda, ...
                basisobj, Kmat, basiscell, inactive)
            n      = length(argvals);
            nbasis = size(cveci,1);
            Wfd    = fd(cveci, basisobj);
            [f,     basiscell]  =   monfn(argvals, Wfd, basiscell);
            [Dyhat, basiscell]  = mongrad(argvals, Wfd, basiscell);
            if ~isempty(zmat)
                xmat = [zmat,f];
            else
                xmat = [ones(n,1), f];
            end
            ncov   = size(xmat,2);
            matwrd = all(size(wtvec) == n);
            if matwrd
                betmi = (xmat'*wtvec*xmat)\(xmat'*wtvec*yi);
            else
                wtroot = sqrt(wtvec);
                wtrtmt = wtroot*ones(1,ncov);
                yroot  = yi.*wtroot;
                xroot  = xmat.*wtrtmt;
                %  compute regression coefs.
                betmi = xroot\yroot;
            end
            %  update fitted values
            yhat  = xmat*betmi;
            %  update residuals and function values
            res     = yi - yhat;
            if matwrd
                Fstri.f = (res'*wtvec*res)./n;
            else
                Fstri.f = mean(res.^2.*wtvec);
            end
            grad    = zeros(nbasis,1);
            Dxmat   = zeros([n,ncov,nbasis]);
            Dxmat(:,ncov,:) = Dyhat;
            for j=1:nbasis
                Dxmatj = squeeze(Dxmat(:,:,j));
                if matwrd
                    yDx = yi'*wtvec*Dxmatj*betmi;
                    xDx = xmat'*wtvec*Dxmatj;
                else
                    Dxroot  = squeeze(Dxmat(:,:,j)).*wtrtmt;
                    yDx     = yroot'*Dxroot*betmi;
                    xDx     = xroot'*Dxroot;
                end
                grad(j) = betmi'*(xDx+xDx')*betmi - 2*yDx;
            end
            Fstri.grad = grad/n;
            if lambda > 0
                Fstri.grad = Fstri.grad +     2 .* Kmat * cveci;
                Fstri.f    = Fstri.f    + cveci' * Kmat * cveci;
            end
            if ~isempty(inactive), Fstri.grad(inactive) = 0; end
            Fstri.norm = sqrt(sum(Fstri.grad.^2));   %  gradient norm
            
            %  ----------------------------------------------------------------
            
            function hessmat = hesscal(betmi, Dyhat, wtvec, lambda, Kmat, inactive)
                nbet = length(betmi);
                [n, nbasis] = size(Dyhat);
                temp = betmi(nbet).*Dyhat;
                matwrd = all(size(wtvec) == n);
                if matwrd
                    hessmat = 2.*temp'*wtvec*temp./n;
                else
                    wtroot = sqrt(wtvec);
                    temp = temp.*(wtroot*ones(1,nbasis));
                    hessmat = 2.*temp'*temp./n;
                end
                %  adjust for penalty
                if lambda > 0, hessmat = hessmat + 2.*Kmat; end
                %  adjust for inactive coefficients
                if ~isempty(inactive)
                    hessmat(inactive,:    ) = 0;
                    hessmat(:       ,inactive) = 0;
                    hessmat(inactive,inactive) = eye(ninactive);
                end
                

