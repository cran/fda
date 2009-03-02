function [regfd, Wfd, shift, Fstr, iternum] = ...
              registerfd(y0fd, yfd, Wfd0Par, periodic,   ...
                         crit, conv, iterlim, dbglev)
%REGISTERFD registers a set of curves YFD to a target function Y0FD.
%  Arguments are:
%  Y0FD    ... Functional data object for target function.  It must
%                contain a single curve, but this single curve
%                can be multivariate.
%  YFD     ... Functional data object for functions to be registered
%  WFD0PAR ... Functional parameter object for function W defining warping fns
%              Its coefficients are the starting values used in the
%                iterative computation of the final warping fns.
%              NB:  The first coefficient is NOT used.
%                For both B-spline and Fourier bases, this first 
%                coefficient determines the constant term in the expansion,
%                and, since a register function is normalized, this term
%                is, in effect, eliminated or has no influence on the 
%                result.  This first position is used, however, to 
%                contain the shift parameter in case the data are 
%                treated as periodic.  At the end of the calculations,
%                the shift parameter is returned separately.
%  PERIODIC... If one, curves are periodic and a shift parameter is fit.
%              Initial value for shift parameter is taken to be 0.  
%              The periodic option should ONLY be used with a Fourier 
%              basis for the target function Y0FD, the functions to be
%              registered, YFD, and the functions WFD0 defining the 
%              time-warping functions.
%  CRIT    ... if 1 least squares, if 2 log eigenvalue ratio.  Default is 1.
%                Default:  0
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  DBGLEV  ... level of output of computation history
%  Returns:
%  REGFD  ...  A functional data object for the registered curves
%  WFD    ...  A Functional data object for function W defining 
%                         warping fns
%  SHIFT  ...  Shift parameter value if curves are periodic

%  last modified 20 July 2006

%  set default arguments

if nargin <  8, dbglev  =  1;    end
if nargin <  7, iterlim = 20;    end
if nargin <  6, conv    = 1e-2;  end
if nargin <  5, crit    =  2;    end
if nargin <  4, periodic=  0;    end
if nargin <  3
   basis  = getbasis(y0fd);
   rangex = getbasisrange(basis);
   if periodic
      wbasis = create_fourier_basis(rangex, 5);
   else
      wbasis = create_bspline_basis(rangex, 5, 4);
   end
   Wfd0Par   = fdPar(wbasis);
end

%  check target function(s)

y0dim  = size(getcoef(y0fd));
ndimy0 = length(y0dim);
if y0dim(ndimy0) == 1, ndimy0 = ndimy0 - 1; end;
if ndimy0 <= 1, nvar = 1; else nvar = y0dim(ndimy0); end
if ndimy0 == 3
   if y0dim(2) > 1
      error('Y0FD is not a single function.');
   end
end
basis0  = getbasis(y0fd);
rangex0 = getbasisrange(basis0);

%  check functions to be registered

ydim = size(getcoef(yfd));
ncurve = ydim(2);
ndimy = length(ydim);
if ydim(ndimy) == 1, ndimy = ndimy - 1; end;
if ndimy == 1
   if ndimy0 > 1, error('YFD is not compatible with Y0FD.'); end
end
if ndimy == 2
   if ndimy0 > 2, error('YFD is not compatible with Y0FD.'); end
end
if ndimy == 3
   if ndimy0 < 3, error('YFD is not compatible with Y0FD.'); end
end
if ndimy > 3, error('YFD is more than 3-dimensional.'); end
ybasis  = getbasis(yfd);
ynbasis = getnbasis(ybasis);
if periodic
   if ~strcmp(getbasistype(ybasis), 'fourier')
      error('PERIODIC is true, basis not fourier type.');
   end
end

%  check Wfd0Par object

if ~isa_fdPar(Wfd0Par) 
    if isa_fd(Wfd0Par) || isa_basis(Wfd0Par)
        Wfd0Par = fdPar(Wfd0Par);
    else
        error(['WFD0PAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up WFD0

Wfd0 = getfd(Wfd0Par);

%  set up LFDOBJ

Lfdobj = getLfd(Wfd0Par);
Lfdobj = int2Lfd(Lfdobj);

%  set up LAMBDA

lambda = getlambda(Wfd0Par);

%  check functions W defining warping functions

wcoef  = getcoef(Wfd0);
wbasis = getbasis(Wfd0);
nbasis = getnbasis(wbasis);
norder = nbasis - length(getbasispar(wbasis));
rangex = getbasisrange(wbasis);
wdim   = size(wcoef);
ncoef  = wdim(1);
if wdim(2) == 1
    wcoef = repmat(wcoef,1,ncurve);
    Wfd0  = putcoef(Wfd0, wcoef);
    wdim  = size(wcoef);
end
ndimw  = length(wdim);
if ndimw == 2
   if wdim(2) ~= ncurve
      error('WFD and YFD do not have the same dimensions.');
   end
end
if ndimw > 2
   error('WFD is not univariate.');
end

%  set up a fine mesh of argument values

NFINEMIN = 101;
nfine = 10*ynbasis + 1;
if nfine < NFINEMIN, nfine = NFINEMIN;  end

xlo   = rangex(1);
xhi   = rangex(2);
width = xhi - xlo;
xfine = linspace(xlo, xhi, nfine)';

%  evaluate target curve at fine mesh of values

xlo0   = rangex0(1);
xhi0   = rangex0(2);
xfine0 = linspace(xlo0, xhi0, nfine)';
y0fine = squeeze(eval_fd(y0fd, xfine0));

%  set up indices of coefficients that will be modified in ACTIVE

wcoef1   = wcoef(1,:);
if periodic
   active   = 1:nbasis;
   wcoef(1) = 0;
else
   active = 2:nbasis;
end
  
%  initialize matrix Kmat defining penalty term

if lambda > 0
   Kmat = eval_penalty(wbasis, Lfdobj);
   ind  = 2:ncoef; 
   Kmat = lambda.*Kmat(ind,ind);
else
   Kmat = [];
end

%  set up limits on coefficient sizes

climit = 50.*[-ones(1,ncoef); ones(1,ncoef)];

%  set up cell for storing basis function values
  
JMAX = 15;
basiscell = cell(1,JMAX);

yregcoef = getcoef(yfd);

%  get penalty matrix for curves to be registered

penmat  = eval_penalty(ybasis);
penmat  = penmat + 1e-10 .* max(max(penmat)) .* eye(ynbasis);
penmat  = sparse(penmat);

%  iterate through the curves

for icurve = 1:ncurve
  if dbglev >= 1 && ncurve > 1
      fprintf(['\n\n-------  Curve ',num2str(icurve),'  --------\n'])
  end;
  if ncurve == 1
    yfdi = yfd;
    Wfdi = Wfd0;
    cvec = wcoef;
  else
    Wfdi = Wfd0(icurve);
    cvec = wcoef(:,icurve);
    if nvar == 1
      yfdi = yfd(icurve);
    else
      yfdi = yfd(icurve,:);
    end
  end
  
  %  evaluate curve to be registered at fine mesh 
  
  yfine = squeeze(eval_fd(xfine, yfdi));
    
  %  evaluate objective function for starting coefficients
  
  %  first evaluate warping function and its derivative at fine mesh
  
  ffine  =   monfn(xfine, Wfdi, basiscell);
  Dffine = mongrad(xfine, Wfdi, basiscell);
  fmax   = ffine(nfine);
  Dfmax  = Dffine(nfine,:);
  hfine  = xlo + width.*ffine./fmax;
  Dhfine = width.*(fmax.*Dffine - ffine*Dfmax)./fmax^2;
  hfine(1)     = xlo;
  hfine(nfine) = xhi;
  
  %  register curves given current Wfdi
  
  if all(cvec == 0)
      yregfdi = yfdi;
  else
      yregfdi = regyfn(xfine, yfine, hfine, yfdi, Wfdi, ...
                       penmat, periodic);
  end
   
  %  compute initial criterion value and gradient
                
  Fstr = regfngrad(xfine, y0fine, Dhfine, yregfdi, Wfdi, ...
                   Kmat, periodic, crit);
  
  %  compute the initial expected Hessian
  
  if crit == 2 
     D2hwrtc = monhess(xfine, Wfdi, basiscell);
     D2fmax  = D2hwrtc(nfine,:);
     m = 0;
     for j=1:nbasis
         khi = min(nbasis,j+norder-1);
         for k=j:khi
             m = m + 1;
             D2hwrtc(:,m) = width.*(2.*ffine.*Dfmax(j).*Dfmax(k)      ...
                - fmax.*(Dffine(:,j).*Dfmax(k) + Dffine(:,k).*Dfmax(j)) ...
                + fmax.^2.*D2hwrtc(:,m) - ffine.*fmax.*D2fmax(m))./fmax.^3;
        end
     end
  else
     D2hwrtc = [];
  end
    
  hessmat = reghess(xfine, y0fine, Dhfine, D2hwrtc, yregfdi, Kmat, ...
                    nbasis, norder, periodic, crit);

  %  evaluate the initial update vector for correcting the initial cvec

  deltac = linesearch(Fstr, hessmat, dbglev);

  %  initialize iteration status arrays

  iternum = 0;
  status = [iternum, Fstr.f, Fstr.norm];
  if dbglev >= 1
    fprintf('\nIter.    Criterion   Grad Length\n')
    fprintf('%3.f %12.4f %10.4f\n', status);
  end
  iterhist = zeros(iterlim+1,length(status));
  iterhist(1,:)  = status;
  if iterlim == 0, break;  end

  %  -------  Begin main iterations  -----------

  MAXSTEPITER = 5;
  MAXSTEP = 100;
  trial   = 1;
  reset   = 0;
  linemat = zeros(3,5);
  cvecold = cvec;
  Foldstr = Fstr;
  dbgwrd  = dbglev >= 2;
  
  %  ---------------  beginning of optimization loop  -----------
  
  for iter = 1:iterlim
      iternum = iternum + 1;
      %  set logical parameters
      dblwrd = [0,0]; limwrd = [0,0]; ind = 0;  ips = 0;
      %  compute slope
      linemat(2,1) = sum(deltac.*Foldstr.grad);
      %  normalize search direction vector
      sdg          = sqrt(sum(deltac.^2));
      deltac       = deltac./sdg;
      linemat(2,1) = linemat(2,1)/sdg;
      % initialize line search vectors
      linemat(:,1:4) = [0; linemat(2,1); Fstr.f]*ones(1,4);
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
      if linemat(2,1) >= -min([1e-3,conv]);
          if dbglev >= 2, disp('Initial slope too small'); end
          wcoef(:,icurve)    = cvec;
          status = [iternum, Fstr.f, Fstr.norm];
          if dbglev >= 1
              fprintf('%3.f %12.4f %10.4f\n', status);
          end
          break;
      end
      %  first step set to trial
      linemat(1,5)  = trial;
      %  ------------  begin line search iteration loop  ----------
      cvecnew = cvec;
      Wfdnewi  = Wfdi;
      for stepiter = 1:MAXSTEPITER
        %  check the step size and modify if limits exceeded
        [linemat(1,5), ind, limwrd] = ...
           stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
                   climit, active, dbgwrd);
        if ind == 1, break; end % break of limit hit twice in a row
        if linemat(1,5) <= 1e-7
           %  Current step size too small ... terminate
           if dbglev >= 2
             fprintf('Stepsize too small: %15.7f\n', linemat(1,5));
           end
           break;
        end
        %  update parameter vector
        cvecnew = cvec + linemat(1,5).*deltac;
        %  compute new function value and gradient
        Wfdnewi = putcoef(Wfdi, cvecnew);
        %  first evaluate warping function and its derivative at fine mesh
        cvectmp = cvecnew;
        cvectmp(1) = 0;
        Wfdtmpi = putcoef(Wfdi, cvectmp);
        ffine  =   monfn(xfine, Wfdtmpi, basiscell);
        Dffine = mongrad(xfine, Wfdtmpi, basiscell);
        fmax   = ffine(nfine);
        Dfmax  = Dffine(nfine,:);
        hfine  = xlo + width.*ffine./fmax;
        Dhfine = width.*(fmax.*Dffine - ffine*Dfmax)./fmax^2;
        hfine(1)     = xlo;
        hfine(nfine) = xhi;
        %  register curves given current Wfdi
        yregfdi = regyfn(xfine, yfine, hfine, yfdi, Wfdnewi, penmat, periodic);
        Fstr = regfngrad(xfine, y0fine, Dhfine, yregfdi, Wfdnewi, Kmat, periodic, crit);
        linemat(3,5) = Fstr.f;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*Fstr.grad);
        if dbglev >= 2
          fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                  [stepiter, linemat(:,5)']);
        end
        %  compute next line search step, also testing for convergence
        [linemat, ips, ind, dblwrd] = ...
              stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbglev);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 || ind == 5, break; end
     end
     %  ------------  end line search iteration loop  ----------
     cvec   = cvecnew;
     Wfdi   = Wfdnewi;
     %  test for function value made worse
     if Fstr.f > Foldstr.f
        %  Function value worse ... warn and terminate
        if dbglev >= 2
          fprintf('Criterion increased, terminating iterations.\n');
          fprintf('%10.4f %10.4f\n',[Foldstr.f, Fstr.f]);
        end
        %  reset parameters and fit
        cvec   = cvecold;
        Wfdi   = putcoef(Wfdi, cvecold);
        Fstr   = Foldstr;
        deltac = -Fstr.grad;
        if dbglev > 2
          for i = 1:nbasis, fprintf('%10.4f%', cvec(i)); end
          fprintf('\n');
        end
        if reset == 1
           %  This is the second time in a row that this
           %     has happened ...  quit
           if dbglev >= 2
             fprintf('Reset twice, terminating.\n');
           end
           break;
        else
           reset = 1;
        end
     else
        %  function value has not increased,  check for convergence
        if abs(Foldstr.f-Fstr.f) < conv
           wcoef(:,icurve)    = cvec;
           status = [iternum, Fstr.f, Fstr.norm];
           if dbglev >= 1
              fprintf('%3.f %12.4f %10.4f\n', status);
           end
           break;
        end
        %  update old parameter vectors and fit structure
        cvecold = cvec;
        Foldstr = Fstr;
        %  update the expected Hessian
        if (iter < iterlim)  
            if crit == 2 
               cvectmp = cvec;
               cvectmp(1) = 0;
               Wfdtmpi = putcoef(Wfdi, cvectmp);
               D2hwrtc = monhess(xfine, Wfdtmpi, basiscell);
               D2fmax  = D2hwrtc(nfine,:);
               %  normalize 2nd derivative
               m = 0;
               for j=1:nbasis
                  for k=j:min(nbasis,j+norder-1)
                     m = m + 1;
                     D2hwrtc(:,m) = width.*(2.*ffine.*Dfmax(j).*Dfmax(k)      ...
                    - fmax.*(Dffine(:,j).*Dfmax(k) + Dffine(:,k).*Dfmax(j)) ...
                    + fmax.^2.*D2hwrtc(:,m) - ffine.*fmax.*D2fmax(m))./fmax.^3;
                  end
               end
            else
               D2hwrtc = [];
            end
            hessmat = reghess(xfine, y0fine, Dhfine, D2hwrtc, yregfdi, Kmat, ...
                              nbasis, norder, periodic, crit);
            %  update the line search direction vector
            deltac = linesearch(Fstr, hessmat, dbglev);
            reset = 0;
        end
     end
     status = [iternum, Fstr.f, Fstr.norm];
     iterhist(iter+1,:) = status;
     if dbglev >= 1
       fprintf('%3.f %12.4f %10.4f\n', status);
     end
  end
  %  ---------------  end of optimization loop  -----------
  wcoef(:,icurve) = cvec;
  if nvar == 1
     yregcoef(:,icurve)   = getcoef(yregfdi);
  else
     yregcoef(:,icurve,:) = getcoef(yregfdi);
  end
end

%  --------------------   end of variable loop  -----------

%  create functional data objects for the registered curves

regfdnames    = getnames(yfd);
regfdnames{3} = ['Registered ',regfdnames{3}];
ybasis        = getbasis(yfd);
regfd         = fd(yregcoef, ybasis, regfdnames);

%  create functional data objects for the warping functions

if periodic
  shift      = wcoef(1,:)';
  wcoef(1,:) = wcoef1;
else
   shift = zeros(ncurve,1);
end
Wfd      = fd(wcoef, wbasis);

%  ----------------------------------------------------------------

function Fstr = regfngrad(xfine, y0fine, Dhfine, yregfd, Wfd, ...
                          Kmat, periodic, crit)
  
  nvar  = size(y0fine,2);
  nfine = length(xfine);
  cvec  = getcoef(Wfd);
  ncvec = length(cvec);
  onecoef = ones(1,ncvec);
  
  if periodic
     Dhfine(:,1) = 1;
  else
     Dhfine(:,1) = 0;  
  end
  yregmat  = squeeze(eval_fd(yregfd, xfine));
  Dyregmat = squeeze(eval_fd(yregfd, xfine, 1));
  
  %  loop through variables computing function and gradient values
  
  Fval = 0;
  gvec = zeros(ncvec,1);
  for ivar = 1:nvar
    y0ivar  =   y0fine(:,ivar);
    ywrthi  =  yregmat(:,ivar);
    Dywrthi = Dyregmat(:,ivar);
    aa = mean(y0ivar.^2);
    cc = mean(ywrthi.^2);
    bb = mean(y0ivar.*ywrthi);
    ff = aa - cc;
    dd = sqrt(ff^2 + 4*bb^2);
    Dywrtc  = (Dywrthi * onecoef).*Dhfine;
    if crit == 1
      res  = y0ivar - ywrthi;
      Fval = Fval + aa - 2*bb + cc;
      gvec = gvec - 2.*Dywrtc'*res./nfine;
    else
      Fval = Fval + aa + cc - dd;
      Dbb  =    Dywrtc'*y0ivar./nfine;
      Dcc  = 2.*Dywrtc'*ywrthi./nfine;
      Ddd  = (4.*bb.*Dbb - ff.*Dcc)./dd;
      gvec = gvec + Dcc - Ddd;
    end
  end
  if ~isempty(Kmat)
     ind   = 2:ncvec;
     ctemp = cvec(ind,1);
     Kctmp = Kmat*ctemp;
     Fval  = Fval + ctemp'*Kctmp;
     gvec(ind) = gvec(ind) + 2.*Kctmp;
  end
  
%  set up F structure containing function value and gradient

  Fstr.f    = Fval;
  Fstr.grad = gvec;
  %  do not modify initial coefficient for B-spline and Fourier bases
  if ~periodic,  Fstr.grad(1) = 0;  end
  Fstr.norm = sqrt(sum(Fstr.grad.^2));

%  ---------------------------------------------------------------

function hessmat = reghess(xfine, y0fine, Dhfine, D2hwrtc, yregfd, ...
                           Kmat, nbasis, norder, periodic, crit)
                        
  nvar      = size(y0fine,2);
  nfine     = length(xfine);
  onecoef   = ones(1,nbasis);
  nbaspr    = nbasis*norder - norder*(norder-1)/2;
  onenbaspr = ones(1,nbaspr);
   
  if periodic
     Dhfine(:,1) = 1;
  else
     Dhfine(:,1) = 0;  
  end
  yregmat  = squeeze(eval_fd(yregfd, xfine));
  Dyregmat = squeeze(eval_fd(yregfd, xfine, 1));
  
  if crit == 2 
     D2yregmat = squeeze(eval_fd(yregfd, xfine, 2));
     %  set values corresponding to first and inactive basis function
     if periodic
        D2hwrtc(:,1) = 0;
        for j=2:norder,  D2hwrtc(:,j) = Dhfine(:,j); end
     else
        D2hwrtc(:,1) = 1;
        for j=2:norder,  D2hwrtc(:,j) = 0;           end
     end   
  end

  hessmat = zeros(nbasis);
  for ivar = 1:nvar
    y0i    =   y0fine(:,ivar);
    yri    =  yregmat(:,ivar);
    D1yri  = Dyregmat(:,ivar);
    Dywrtc = ((D1yri * onecoef).*Dhfine);
    if crit == 1
      hessmat = hessmat + 2.*(Dywrtc' * Dywrtc)./nfine;
    else
      D2yri = D2yregmat(:,ivar);
      aa      = mean(y0i.^2);
      cc      = mean(yri.^2);
      bb      = mean(y0i.*yri);
      ff      = aa - cc;
      dd      = sqrt(ff.*ff + 4*bb*bb);
      Dbb     =    Dywrtc'*y0i./nfine;
      Dcc     = 2.*Dywrtc'*yri./nfine;
      wgtbb1  = (y0i.*D2yri           )*onecoef;
      wgtcc1  = (yri.*D2yri + D1yri.^2)*onecoef;
      D2bb    =    (Dhfine.*wgtbb1)'*Dhfine./nfine;
      D2cc    = 2.*(Dhfine.*wgtcc1)'*Dhfine./nfine;
      wgtbb2  = (y0i.*D1yri)*onenbaspr;
      wgtcc2  = (yri.*D1yri)*onenbaspr;
      y0DyD2h =    mean(wgtbb2.*D2hwrtc)';
      yrDyD2h = 2.*mean(wgtcc2.*D2hwrtc)';
      m = 0;
      for j=1:nbasis
          for k=j:min(nbasis,j+norder-1)
              m = m + 1;
              D2bb(j,k) = D2bb(j,k) + y0DyD2h(m);
              D2bb(k,j) = D2bb(j,k);
              D2cc(j,k) = D2cc(j,k) + yrDyD2h(m);
              D2cc(k,j) = D2cc(j,k);
          end
      end
      DbbDcc = Dbb*Dcc';
      term1  = - 4.*bb.*D2bb./dd;
      term2  = (1 + ff/dd).*D2cc;
      term3  = -((2+(4*bb/dd))*(2-(4*bb/dd)).*Dbb*Dbb' + ...
                 (4*ff*bb/dd^2).*(DbbDcc + DbbDcc')    + ...
                 (1+(ff/dd))*(1-(ff/dd)).*Dcc*Dcc')./dd;
      hessmat = hessmat + term1 + term2 + term3; 
    end
  end
  if ~isempty(Kmat)
     ind = 2:nbasis;
     hessmat(ind,ind) = hessmat(ind,ind) + 2.*Kmat;
  end
  if ~periodic
     hessmat(1,:) = 0;
     hessmat(:,1) = 0;
     hessmat(1,1) = 1;
  end
  %hesseig  = sort(eig(hessmat));
  %disp( hesseig([1, 2, 3])')
  %disp( hesseig([nbasis-2, nbasis-1, nbasis])')
  
%  ----------------------------------------------------------------

function yregfd = regyfn(xfine, yfine, hfine, yfd, Wfd, penmat, ...
                         periodic)

%  get shift value for the periodic case from Wfd

coef   = getcoef(Wfd);
shift  = coef(1);  
coef(1)= 0;

%  if all coefficients are zero, no transformation  needed

if all(coef == 0)
   if periodic
      if shift == 0
         yregfd = yfd;
         return;
      end
   else
      yregfd = yfd;
      return;
   end
end

%  Estimate inverse of warping function at fine mesh of values  
%  28 dec 000
%  It makes no real difference which 
%     interpolation method is used here.
%  Linear is faster and sure to be monotone.
%  Using WARPSMTH added nothing useful, and was abandoned.

nfine       = length(xfine);
hinv        = safeinterp(hfine, xfine, xfine);
hinv(1)     = xfine(1);
hinv(nfine) = xfine(nfine);

%  carry out shift if period and shift ~= 0

if periodic && shift ~= 0
   yfine = shifty(xfine, yfine, shift);
end

%  smooth relation between Y and HINV
%  this is the same code as in PROJECT_BASIS, but avoids
%  recomputing the penalty matrix

basis    = getbasis(yfd);
if any(hinv < xfine(1))
    save hfine
    save xfine
    save hinv
    error(['HINV values too small by ',num2str(xfine(1)-min(hinv))]);
end
if any(hinv > xfine(nfine))
    error(['HINV values too large by ',num2str(max(hinv)-xfine(nfine))]);
end
basismat = getbasismatrix(hinv, basis);
Bmat     = basismat' * basismat;
lambda1  = (0.0001 .* sum(diag(Bmat)))./sum(diag(penmat));
Cmat     = Bmat + lambda1 .* penmat;
Dmat     = basismat' * yfine;
ycoef    = symsolve(Cmat,Dmat);

%  set up FD object for registered function

yregfd   = fd(ycoef, basis);

%  ----------------------------------------------------------------

function [deltac, cosangle] = linesearch(Fstr, hessmat, dbglev)
ncoef   = size(hessmat,1);
hesseig = sort(eig(hessmat));
if hesseig(1) == 0
    delta = 1e-5*hesseig(ncoef);
    hessmat = hessmat + delta.*eye(ncoef);
end
if hesseig(1) < 0
    delta = 1e-5*hesseig(ncoef) - hesseig(1);
    hessmat = hessmat + delta.*eye(ncoef);
end
deltac = -symsolve(hessmat, Fstr.grad);
cosangle = -sum(Fstr.grad.*deltac)./sqrt(sum(Fstr.grad.^2)*sum(deltac.^2));
if dbglev >= 2
   fprintf('Cos(angle) = %10.4f\n',cosangle);
end
if cosangle < 1e-7 
   if dbglev >=2, fprintf('angle negative\n'); end
   deltac = -Fstr.grad;
end


