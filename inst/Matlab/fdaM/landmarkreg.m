function [regfd, warpfd, Wfd] = landmarkreg(fdobj, ximarks, x0marks, ...
                                            WfdParobj, monwrd, ylambda)
%  LANDMARKREG ... Register curves using landmarks.
%  Arguments:
%  FDOBJ     ...  a functional data object for curves to be registered
%  XIMARKS   ...  N by NL array of times of interior landmarks for
%                 each observed curve
%  XOMARKS   ...  vector of length NL of times of interior landmarks for
%                 target curve
%  WFDPAROBJ ...  A functional parameter object used to define the
%                 strictly monotone warping function.
%  MONWRD    ...  If 1, warping functions are estimated by monotone smoothing,
%                 otherwise by regular smoothing.  The latter is faster, but
%                 not guaranteed to produce a strictly monotone warping 
%                 function.  If MONWRD is 0 and an error message results 
%                 indicating nonmonotonicity, rerun with MONWRD = 1.
%                 Default:  1
%  YLAMBDA ... smoothing parameter to be used in computing the registered
%                 functions.  For high dimensional bases, local wiggles may be 
%                 found in the registered functions or its derivatives that are
%                 not seen in the unregistered functions.  In this event, this
%                 parameter should be increased until they disappear.             
%  Returns:
%  REGFD  ...  A functional data object for the registered curves
%  WARPFD ...  A functional data object for the warping functions.  
%              This tends only to be useful if monwrd = 0.
%  WFD    ...  A Functional data object for function W defining 
%              warping fns.  This is only produced for monwrd = 1.
%              The warping function values for case i require code 
%              of the form
%                     harg    = monfn(evalarg, Wfd(i));
%                     width   = rng(2) - rng(1);
%                     warpvec = rng(1) + width.*harg./rng(2);
%              where rng is the argument range vector.

%  last modified 15 May 2009 by Jim Ramsay

%  check FDOBJ

if ~isa_fd(fdobj)
    error ('Argument FDOBJ is not a functional data object.');
end

%  extract information from curve functional data object and its basis

coef   = getcoef(fdobj);
coefd  = size(coef);
ndim   = length(coefd);
ncurve = coefd(2);
if ndim > 2
    nvar = coefd(3);
else
    nvar = 1;
end

fdbasisobj = getbasis(fdobj);
fdnbasis   = getnbasis(fdbasisobj);
rangeval   = getbasisrange(fdbasisobj);
fdParobj   = fdPar(fdobj, 2, ylambda);

%  check landmarks

ximarksd = size(ximarks);
if ximarksd(1) ~= ncurve
    error('Number of rows of second argument wrong.');
end
nlandm = ximarksd(2);

if nargin < 3
    x0marks = mean(ximarks);
end
if (length(x0marks) ~= nlandm)
    error(['Number of target landmarks not equal to ', ...
            'number of curve landmarks.']);
end
if size(x0marks,2) == 1, x0marks = x0marks';  end

%  set default argument values

if nargin < 6, ylambda = 1e-10;  end
if nargin < 5, monwrd  = 1;      end

if nargin < 4
    basisobj  = getbasis(fdobj);
    rangex    = getbasisrange(basisobj);
    wnbasis   = length(x0marks) + 2;
    wbasis    = create_bspline_basis(rangex, wnbasis);
    WfdParobj = fdPar(wbasis);
else
    %  check  WfdParobj object
    WfdParobj = fdParcheck(WfdParobj);
end

%  set up WFD0 and WBASIS

Wfd0   = getfd(WfdParobj);
wbasis = getbasis(Wfd0);

%  set up LFDOBJ

WLfdobj = getLfd(WfdParobj);
WLfdobj = int2Lfd(WLfdobj);

%  set up LAMBDA

wlambda = getlambda(WfdParobj);

%  check landmark target values

test1 = any(x0marks(:) <= rangeval(1));
test2 = any(ximarks(:) >= rangeval(2));
if  test1 || test2
    error('Some landmark values are not within the range.');
end

%  set up analysis

n   = min([101,10*fdnbasis]);
x   = linspace(rangeval(1),rangeval(2),n)';
y   = eval_fd(x, fdobj);
yregmat = y;
hfunmat = zeros(n,ncurve);
wlambda = max([wlambda,1e-10]);

xval      = [rangeval(1),x0marks,rangeval(2)]';
nwbasis   = getnbasis(wbasis);
Wfd0      = fd(zeros(nwbasis,1),wbasis);
Wcoef     = zeros(nwbasis,ncurve);
WfdParobj = fdPar(Wfd0, WLfdobj, wlambda);

%  --------------------------------------------------------------------
%                  Iterate through curves to register
%  --------------------------------------------------------------------

for icurve = 1:ncurve
    
    %  set up landmark times for this curve
    fprintf(['Curve ',num2str(icurve),'\n']);
    yval = [rangeval(1),ximarks(icurve,:),rangeval(2)]';
    
    if all(ximarks(icurve,:) == x0marks)
        Wcoef(:,icurve)   = zeros(nwbasis,1);
        hfunmat(:,icurve) = x;
        yregmat(:,icurve) = eval_fd(x, fdobj(icurve));
    else
        %  smooth relation between this curve's values and target's values
        %  warpfd is the functional data object for the warping functions
        %  h is a vector of warping function values corresponding to arguments x
        if monwrd
            %  use monotone smoother
            Wfd      = smooth_morph(xval, yval, WfdParobj);
            h        = monfn(x, Wfd);
            %  normalize h
            b        = (rangeval(2)-rangeval(1))/(h(n)-h(1));
            a        = rangeval(1) - b*h(1);
            h        = a + b.*h;
            h([1,n]) = rangeval;
        else
            %  use regular smoother
            warpfd   = smooth_basis(xval, yval, wbasis);
            %  set up warping function sampling values
            h        = eval_fd(x, warpfd);
            %  normalize h
            b        = (rangeval(2)-rangeval(1))/(h(n)-h(1));
            a        = rangeval(1) - b*h(1);
            h        = a + b.*h;
            h([1,n]) = rangeval;
            %  check for monotonicity because regular smooth may not be monotone
            deltah = diff(h);
            if any(deltah <= 0) 
                error(['Non-increasing warping function estimated for curve',...
                        num2str(icurve),'.\n',...
                        'Try setting MONWRD to 1.']);
            end
        end
        hfunmat(:,icurve) = h;
        
        %  compute h-inverse in order to register curves
        
        if monwrd
            Wfdinv      = smooth_morph(h, x, WfdParobj);
            hinv        = monfn(x, Wfdinv);
            b           = (rangeval(2)-rangeval(1))/(hinv(n)-hinv(1));
            a           = rangeval(1) - b*hinv(1);
            hinv        = a + b.*hinv;
            hinv([1,n]) = rangeval;
        else
            hinvfd      = smooth_basis(h, x, WfdParobj);
            hinv        = eval_fd(x, hinvfd);
            b           = (rangeval(2)-rangeval(1))/(hinv(n)-hinv(1));
            a           = rangeval(1) - b*hinv(1);
            hinv        = a + b.*hinv;
            hinv([1,n]) = rangeval;
            %  check for monotonicity because regular smooth may not be 
            %  monotone
            deltahinv = diff(hinv);
            if any(deltahinv <= 0) 
                error(['Non-increasing inverse warping function ', ...
                       'estimated for curve',num2str(icurve)]);
            end
        end
        
        %  compute registered curves
        
        if ndim == 2
            %  single variable case
            yregfd = smooth_basis(hinv, y(:,icurve), fdParobj);
            yregmat(:,icurve) = eval_fd(x, yregfd);
        end
        if ndim == 3
            %  multiple variable case
            for ivar = 1:nvar
                % evaluate curve as a function of h at sampling points
                ymati  = squeeze(y(:,icurve,ivar));
                yregfd = smooth_basis(hinv, ymati, fdParobj);
                yregmat(:,icurve,ivar) = eval_fd(x, yregfd);
            end
        end
    end
end

%  create functional data objects for the registered curves

fdParobj      = fdPar(fdbasisobj, 2, 1e-10);
regfd         = smooth_basis(x, yregmat, fdParobj);
regfdnames    = getnames(fdobj);
regfdnames{3} = ['Registered ',regfdnames{3}];
regfd         = putnames(regfd, regfdnames);

%  create functional data objects for the warping functions

warpfd         = smooth_basis(x, hfunmat, fdParobj);
warpfdnames    = getnames(fdobj);
warpfdnames{3} = ['Warped ',warpfdnames{3}];
warpfd         = putnames(warpfd, warpfdnames);

%  create functional data object for Wfd if monwrd is true

Wfd = fd(Wcoef,wbasis);



