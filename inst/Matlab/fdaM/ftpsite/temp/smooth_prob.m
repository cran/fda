function [Wfd, Fstr, iternum, iterhist] = ...
   smooth_prob(argvals, y, Wfd, Lfdobj, lambda, conv, iterlim, dbglev)
%SMOOTH_PROB smooths the relationship of binomial observations to ARGVALS. 
%  The fitting criterion is penalized log likelihood.

%  Arguments are:
%  ARGVALS ... A vector of argument values
%  Y       ... A vector of data
%  WFD     ... functional data basis object defining initial density
%  LFDOBJ  ... linear differential operator defining roughness penalty
%  LAMBDA  ... smoothing parameter
%  ITERLIM ... iteration limit for scoring iterations
%  CONV    ... convergence criterion
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFD      ... functional data basis object defining final density
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  last modified 20 July 2006

if nargin < 2
   error('Less than two arguments are supplied.');
end

%  ensure that both argvals and y are column vectors

if size(argvals,1) > 1 && size(argvals,2) > 1
    error('Argument X is not a vector.');
end
if size(y,1) > 1 && size(y,2) > 1
    error('Argument Y is not a vector.');
end
argvals = argvals(:);
y = y(:);

%  get basis information

basis  = getbasis(Wfd);
nbasis = getnbasis(basis);
rangex = getbasisrange(basis);
active = 1:nbasis;

%  deal with values out of range

inrng = find(argvals >= rangex(1) && argvals <= rangex(2));
if (length(inrng) ~= length(argvals))
    warning('Some values in X out of range and not used.');
end

argvals     = argvals(inrng);
y     = y(inrng);

%  set some default arguments and constants

if nargin < 9, dbglev  = 1;          end
if nargin < 8, conv    = 1e-4;       end
if nargin < 7, iterlim = 50;         end
if nargin < 6, lambda  = 0;          end
if nargin < 5, Lfdobj  = int2Lfd(2); end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  set up some arrays

climit  = [-50,0;0,400]*ones(2,nbasis);
cvec0   = getcoef(Wfd);
dbgwrd  = dbglev > 1;

%  initialize matrix Kmat defining penalty term

if lambda > 0
  Kmat = lambda.*eval_penalty(basis, Lfdobj);
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl] = loglfnpos(argvals, y, basis, cvec0);

%  compute initial badness of fit measures

Foldstr.f  = -logl;
gvec  = -Dlogl;
if lambda > 0
   gvec = gvec + 2.*(Kmat * cvec0);
   Foldstr.f = Foldstr.f + cvec0' * Kmat * cvec0;
end
Foldstr.norm = sqrt(mean(gvec.^2));

%  compute the initial expected Hessian

hmat = loglhesspos(argvals, y, basis, cvec0);
if lambda > 0
    hmat = hmat + 2.*Kmat;
end

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
        [logl, Dlogl] = loglfnpos(argvals, y, basis, cvecnew);
        Fstr.f  = -logl;
        gvecnew = -Dlogl;
        if lambda > 0
            gvecnew = gvecnew + 2.*Kmat * cvecnew;
            Fstr.f = Fstr.f + cvecnew' * Kmat * cvecnew;
        end
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
            stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbgwrd);
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
    hmat = loglhesspos(argvals, y, basis, cvec);
    if lambda > 0
        hmat = hmat + 2.*Kmat;
    end
    %  evaluate the update vector
    deltac = -hmat\gvec;
    cosangle  = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));
    if cosangle < 0
        if dbglev > 1, disp('cos(angle) negative'); end
        deltac = -gvec;
    end
    Foldstr = Fstr;
end

%  ---------------------------------------------------------------

function [logl, Dlogl] = loglfnpos(argvals, y, basis, cvec)
N = length(argvals);
phimat  = getbasismatrix(argvals, basis);
Wvec    = phimat*cvec;
EWvec   = exp(Wvec);
res     = y - EWvec;
logl    = -mean(res.^2);
Dlogl   = 2.*phimat'*(res.*EWvec)./N;

%  ---------------------------------------------------------------

function D2logl = loglhesspos(argvals, y, basis, cvec)
N = length(argvals);
nbasis  = getnbasis(basis);
phimat  = getbasismatrix(argvals, basis);
Wvec    = phimat*cvec;
EWvec   = exp(Wvec);
res     = y - EWvec;
Dres    = ((res.*EWvec)*ones(1,nbasis)) .* phimat;
D2logl   = 2.*Dres'*Dres./N;

