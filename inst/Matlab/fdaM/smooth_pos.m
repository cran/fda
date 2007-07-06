function [Wfd, Fstr, iternum, iterhist] = ...
    smooth_pos(argvals, y, fdParobj, wtvec, conv, iterlim, dbglev)
%SMOOTH_POS smooths the relationship of Y to ARGVALS by fitting 
%     a positive function of the form
%                   f(t) = exp W(t)
%     where  W  is a function defined over the same range as ARGVALS,
%  W  = log f.
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum w_i[y_i - f(t_i)]^2 +
%                     \lambda * \int [L W]^2 
%  The function W(t) is expanded by the basis in functional data object
%    Wfd.   The coefficients of this expansion are called "coefficients"
%    in the comments, while the b's are called "regression coefficients"

%  Arguments are:
%  ARGVALS ... A vector of argument values
%  Y       ... A vector of data values
%  FDPAROBJ ... A functional parameter or fdPar object.  This object 
%               contains the specifications for the functional data
%               object to be estimated by smoothing the data.  See
%               comment lines in function fdPar for details.
%               The functional data object WFD in FDPAROBJ is used
%               to initialize the optimization process.
%               It's coefficient array has a single column, and these 
%               are the starting values for the iterative minimization 
%               of mean squared error.
%               This argument may also be either a FD object, or a 
%               BASIS object.  In this case, the smoothing parameter 
%               LAMBDA is set to 0.
%  WTVEC   ...  Vector of positive weights for observations
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFD      ... functional data basis object defining final density
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  last modified 16 June 2007

%  set some default arguments and constants

N = length(y);

if nargin < 7, dbglev  = 1;          end
if nargin < 6, iterlim = 20;         end
if nargin < 5, conv    = 1e-4;       end
if nargin < 4, wtvec   = ones(N,1);  end
if nargin < 3
    error('Less than three arguments are supplied.');
end

%  check ARGVALS

if ~strcmp(class(argvals), 'double')
    error('ARGVALS is not of class double.');
end

if size(argvals,1) == 1
    argvals = argvals';
end

[nobs, ncl] = size(argvals);  %  number of observations
if ncl > 1
    error('ARGVALS is not a vector.')
end
if nobs < 3
    error('ARGVALS does not contain at least three values.');
end

%  check Y

if ~strcmp(class(y), 'double')
    error('Y is not of class double.');
end

if size(y,1) == 1
    y = y';
end

[nrw, ncl] = size(y);  %  number of observations

if ncl > 1
    error('Y is not a vector.')
end
if nrw ~= nobs
    error('Y is not the same length as ARGVALS.');
end

%  check fdParobj

if ~isa_fdPar(fdParobj) 
    if isa_fd(fdParobj) || isa_basis(fdParobj)
        fdParobj = fdPar(fdParobj);
    else
        error(['FDPAROBJ is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up LFDOBJ

Lfdobj = getLfd(fdParobj);
Lfdobj = int2Lfd(Lfdobj);

%  set up LAMBDA

lambda = getlambda(fdParobj);

%  set up BASIS

Wfd      = getfd(fdParobj);
basisobj = getbasis(Wfd);
nbasis   = getnbasis(basisobj);
rangex   = getbasisrange(basisobj);
active   = 1:nbasis;

%  deal with values out of range

inrng = find(argvals >= rangex(1) & argvals <= rangex(2));
if (length(inrng) ~= length(argvals))
    warning('Wid:range', ...
        'Some values in X out of range and not used.');
end

argvals = argvals(inrng);
y       = y(inrng);

%  set up some arrays

climit  = [-50,0;0,400]*ones(2,nbasis);
cvec0   = getcoef(Wfd);
dbgwrd  = dbglev > 1;

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda.*eval_penalty(basisobj, Lfdobj);
else
    Kmat = zeros(nbasis);
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl] = loglfnpos(argvals, y, wtvec, basisobj, cvec0, Kmat);

%  compute initial badness of fit measures

gvec         = -Dlogl;
Foldstr.f    = -logl;
Foldstr.norm = sqrt(mean(gvec.^2));

%  compute the initial expected Hessian

hmat = loglhesspos(argvals, y, wtvec, basisobj, cvec0, Kmat);

%  evaluate the initial update vector for correcting the initial bmat

deltac   = -hmat\gvec;

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Foldstr.f, -logl, Foldstr.norm];
fprintf('\nIteration  Criterion  Neg. Log L  Grad. Norm\n')
fprintf('\n%5.f     %10.4f %10.4f %10.4f\n', status);
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:)  = status;
if iterlim == 0
    Fstr = Foldstr;
    iterhist = iterhist(1,:);
    return;
end

%  -------  Begin iterations  -----------

STEPMAX = 5;
MAXSTEP = 400;
trial   = 1;
cvec    = cvec0;
linemat = zeros(3,5);

for iter = 1:iterlim
    iternum = iternum + 1;
    %  take optimal stepsize
    dblwrd = [0,0]; limwrd = [0,0]; ind = 0;
    %  compute slope
    Fstr = Foldstr;
    linemat(2,1) = sum(deltac.*gvec);
    %  normalize search direction vector
    sdg     = sqrt(sum(deltac.^2));
    deltac  = deltac./sdg;
    linemat(2,1) = linemat(2,1)/sdg;
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        disp('Initial slope nonnegative.')
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-5;
        if dbglev>1, disp('Initial slope too small'); end
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    linemat(1,1:4) = 0;
    linemat(2,1:4) = linemat(2,1);
    linemat(3,1:4) = Foldstr.f;
    stepiter  = 0;
    if dbglev > 1
        fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
            [stepiter, linemat(:,1)']);
    end
    ips = 0;
    %  first step set to trial
    linemat(1,5) = trial;
    %  Main iteration loop for linesrch
    for stepiter = 1:STEPMAX
        %  ensure that step does not go beyond limits on parameters
        %  check the step size
        [linemat(1,5),ind,limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
            climit, active, dbgwrd);
        if linemat(1,5) <= 1e-9
            %  Current step size too small ... terminate
            Fstr    = Foldstr;
            cvecnew = cvec;
            gvecnew = gvec;
            if dbglev > 1
                fprintf('Stepsize too small:  %10.4f\n', linemat(1,5));
            end
            break;
        end
        cvecnew = cvec + linemat(1,5).*deltac;
        %  compute new function value and gradient
        [logl, Dlogl] = ...
            loglfnpos(argvals, y, wtvec, basisobj, cvecnew, Kmat);
        gvecnew   = -Dlogl;
        Fstr.f    = -logl;
        Fstr.norm = sqrt(mean(gvecnew.^2));
        linemat(3,5) = Fstr.f;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*gvecnew);
        if dbglev > 1
            fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
                [stepiter, linemat(:,5)']);
        end
        %  compute next step
        [linemat,ips,ind,dblwrd] = ...
                                stepit(linemat, ips, dblwrd, MAXSTEP);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 || ind == 5, break; end
        %  end iteration loop
    end
    
    cvec = cvecnew;
    gvec = gvecnew;
    Wfd  = putcoef(Wfd, cvec);
    status = [iternum, Fstr.f, -logl, Fstr.norm];
    iterhist(iter+1,:) = status;
    fprintf('%5.f     %10.4f %10.4f %10.4f\n', status);
    %  test for convergence
    if abs(Fstr.f-Foldstr.f) < conv
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    if Fstr.f >= Foldstr.f, break; end
    %  compute the Hessian
    hmat = loglhesspos(argvals, y, wtvec, basisobj, cvec, Kmat);
    %  evaluate the update vector
    deltac   = -hmat\gvec;
    cosangle = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));
    if cosangle < 0
        if dbglev > 1, disp('cos(angle) negative'); end
        deltac = -gvec;
    end
    Foldstr = Fstr;
end

%  ---------------------------------------------------------------

function [logl, Dlogl] = ...
    loglfnpos(argvals, y, wtvec, basisobj, cvec, Kmat)
N       = length(argvals);
phimat  = getbasismatrix(argvals, basisobj);
Wvec    = phimat*cvec;
EWvec   = exp(Wvec);
res     = y - EWvec;
logl    = -mean(wtvec.*res.^2) - cvec'*Kmat*cvec;
Dlogl   = 2.*phimat'*(wtvec.*res.*EWvec)./N - 2.*Kmat*cvec;

%  ---------------------------------------------------------------

function D2logl = loglhesspos(argvals, y, wtvec, basisobj, cvec, Kmat)
N = length(argvals);
nbasis  = getnbasis(basisobj);
phimat  = getbasismatrix(argvals, basisobj);
Wvec    = phimat*cvec;
EWvec   = exp(Wvec);
res     = y - EWvec;
Dres    = ((res.*EWvec)*ones(1,nbasis)) .* phimat;
D2logl  = 2.*Dres'*((wtvec*ones(1,nbasis)).*Dres)./N + 2.*Kmat;

