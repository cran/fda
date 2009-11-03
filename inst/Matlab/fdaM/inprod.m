function  [inprodmat, iter] = ...
                     inprod(fdobj1, fdobj2, Lfdobj1, Lfdobj2, rng, wtfd)
%  INPROD   Computes matrix of inner products of functions.
%  These functions can be:
%  (1) regular functional data objects
%  (2) bivariate functional data objects, in which case they are
%      evaluated along the diagonal
%  (3) basis objects, in which case they are converted to functional
%      data objects with the identity matrix as the coefficient matrix.
%  (4) cell objects, in which each cell (normally there are two) contains
%      a functional data object, and the product of these is used as one
%      of the factors in the inner product.
%
%    If both functions have the same B-spline basis, both Lfdobj's are
%       numeric, there is no wtfd argument, and the ranges are the same,
%       the inner products are exact, and computed by inprod_bspline.
%
%    Otherwise, by numerical integration using Romberg integration 
%       with the trapezoidal rule.

%  Arguments:
%  FDOBJ1 and FDOBJ2 ...  these may be either functional data or basis function
%                    objects.  In the latter case, a functional data object
%                    is created from a basis function object by using the
%                    identity matrix as the coefficient matrix.
%                    Both functional data objects must be univariate.
%                    If  inner products for multivariate objects are needed,
%                    use a loop and call inprod(fdobj1(i),fdobj2(i)).
%  LFDOBJ1 and LFDOBJ2 ...  linear differential operators for inner product for
%                    FD1 and FD2, respectively.  Default to 0.  
%  RNG  ...  Limits of integration
%  WTFD ...  A functional data object defining a weight
%  Return:
%  A NREP1 by NREP2 matrix INPRODMAT of inner products for each possible pair
%  of functions.

%  last modified 2 December 2006

%  set up default values of arguments

if nargin < 6, wtfd = 0;               end
if nargin < 4, Lfdobj2 = int2Lfd(0);   end
if nargin < 3, Lfdobj1 = int2Lfd(0);   end

%  check LFDOBJ1 and LFDOBJ2

Lfdobj1 = int2Lfd(Lfdobj1);
Lfdobj2 = int2Lfd(Lfdobj2);
  
%  set constants determining Richardson extrapolation
  
EPS  = 1e-4;  %  convergence criterion
JMAX = 15;    %  maximum number of iterations
JMIN =  5;    %  minimum number of iterations

%  check    WTFD

if  isa_fd(wtfd)
    coefw = getcoef(wtfd);
    coefd = size(coefw);
    if coefd(2) > 1
        error('Argument WTFD is not a single function');
    end
end
  
%  check FDOBJ1 for being a functional data object or
%  a bivariate functional data object or a basis or a cell
%  containg functional data objects.
%  If a basis, convert to a functional data object having
%  the identity matrix as the coefficient matrix.

fdclass = 1;
if  isa_fd(fdobj1) || isa_bifd(fdobj1)
    coef1  = getcoef(fdobj1);
elseif isa_basis(fdobj1)
    coef1  = eye(getnbasis(fdobj1));
    temp1  = fd(coef1, fdobj1);
    fdobj1 = temp1;
elseif iscell(fdobj1)
    % no action needed
else
    fdclass = 0;
end

%  Check FDOBJ2 for being a functional data object or
%  a bivariate functional data object or a basis.
%  If a basis, convert to a functional data object having
%  the identity matrix as the coefficient matrix.

if  isa_fd(fdobj2) || isa_bifd(fdobj2)
    coef2 = getcoef(fdobj2);
elseif isa_basis(fdobj2)
    coef2  = eye(getnbasis(fdobj2));
    temp2  = fd(coef2, fdobj2);
    fdobj2 = temp2;
elseif iscell(fdobj2)
    %  no action needed
else
    fdclass = 0;
end

% if either FDOBJ1 or FDOBJ2 fails these tests, error message

if ~fdclass
    error (['The two first arguments are not', ...
            ' functional data or basis objects.']);
end

%  check components if FDOBJ1 is cell object

if iscell(fdobj1)
    if length(fdobj1) > 2
        error('More than two cells found for a product of FD objects.');
    end
    fdobj11 = fdobj1{1};
    if ~isa_fd(fdobj11)
        error('A factor in a FD product is not a FD object.');
    end
    fdobj12 = fdobj1{2};
    if ~isa_fd(fdobj12)
        error('A factor in a FD product is not a FD object.');
    end
    %  check coefficient matrices for compatibility
    coef11  = getcoef(fdobj11);
    coef12  = getcoef(fdobj12);
    coefd11 = size(coef11);
    coefd12 = size(coef12);
    if length(coefd11) > 2 || length(coefd12) > 2
        error('A factor in an FD product is not a univariate FD object.');
    end
    if coefd11(2) > 1 && coefd12(2) > 1 && coefd11(2) ~= coefd12(2)
        error('Numbers of replications differ for factors in an FD product.');
    end
    nrep1 = max([coefd11(2),coefd12(2)]);
    %  check ranges for compatibility
    basisobj11 = getbasis(fdobj11);
    basisobj12 = getbasis(fdobj12);
    range11 = getbasisrange(basisobj11);
    range12 = getbasisrange(basisobj12);
    if any(range11 ~= range12)
        error('Ranges differ for factors in an FD product.');
    end
    range1 = range11;
end

%  check components if FDOBJ2 is cell object

if iscell(fdobj2)
    if length(fdobj2) > 2
        error('More than two cells found for a product of FD objects.');
    end
    fdobj21 = fdobj2{1};
    if ~isa_fd(fdobj21)
        error('A factor in a FD product is not a FD object.');
    end
    fdobj22 = fdobj2{2};
    if ~isa_fd(fdobj22)
        error('A factor in a FD product is not a FD object.');
    end
    coef21  = getcoef(fdobj21);
    coef22  = getcoef(fdobj22);
    coefd21 = size(coef21);
    coefd22 = size(coef22);
    if length(coefd21) > 2 || length(coefd22) > 2
        error('A factor in an FD product is not a univariate FD object.');
    end
    if coefd21(2) > 1 && coefd22(2) > 1 && coefd21(2) ~= coefd22(2)
        error('Numbers of replications differ for factors in an FD product.');
    end
    nrep2 = max([coefd21(2),coefd22(2)]);
    %  check ranges for compatibility
    basisobj21 = getbasis(fdobj21);
    basisobj22 = getbasis(fdobj22);
    range21 = getbasisrange(basisobj21);
    range22 = getbasisrange(basisobj22);
    if any(range21 ~= range22)
        error('Ranges differ for factors in an FD product.');
    end
end

%  Determine NREP1 and NREP2, and check for common range.
%  This must be done differently dependent on whether the
%  functional data object is univariate or bivariate.

%  FDOBJ1

if isa_fd(fdobj1)
    coefd1 = size(coef1);
    if length(coefd1) > 2
        error('Functional data object must be univariate');
    end
    nrep1 = coefd1(2);
    basisobj1 = getbasis(fdobj1);
    type1  = getbasistype(basisobj1);
    range1 = getbasisrange(basisobj1);
elseif isa_bifd(fdobj1)
    coefd1 = size(coef1);
    if length(coefd1) > 3
        error('Bivariate functional data object must be univariate');
    end
    if length(coefd1) == 2
        nrep1 = 1;
    else
        nrep1 = coefd1(3);
    end
    sbasisobj1 = getsbasis(fdobj1);
    stype1  = getbasistype(sbasisobj1);
    srange1 = getbasisrange(sbasisobj1);
    tbasisobj1 = gettbasis(fdobj1);
    ttype1  = getbasistype(tbasisobj1);
    trange1 = getbasisrange(tbasisobj1);
    if any(srange1 ~= trange1)
        error('Ranges are not equal for bivariate function FDOBJ1.');
    else
        range1 = srange1;
    end
end

%  FDOBJ2

if isa_fd(fdobj2)
    coefd2 = size(coef2);
    if length(coefd2) > 2
        error('Functional data object must be univariate');
    end
    nrep2 = coefd2(2);
    basisobj2 = getbasis(fdobj2);
    type2  = getbasistype(basisobj2);
elseif isa_bifd(fdobj2)
    coefd2 = size(coef2);
    if length(coefd2) > 3
        error('Bivariate functional data object must be univariate');
    end
    if length(coefd2) == 2
        nrep2 = 1;
    else
        nrep2 = coefd2(3);
    end
    sbasisobj2 = getsbasis(fdobj2);
    stype2  = getbasistype(sbasisobj2);
    srange2 = getbasisrange(sbasisobj2);
    tbasisobj2 = gettbasis(fdobj2);
    ttype2  = getbasistype(tbasisobj2);
    trange2 = getbasisrange(tbasisobj2);
    if any(srange2 ~= trange2)
        error('Ranges are not equal for bivariate function FDOBJ2.');
    end
end

%  set default range

if nargin < 5  
    rng = range1; 
end

%  set iter

iter = 0;

% check range

if rng(1) < range1(1) || rng(2) > range1(2)
    error('Limits of integration are inadmissible.');
end
  
%  Call B-spline version if ...
%  (1) both functional data objects are univariate
%  (2) both bases are B-splines
%  (3) the two bases are identical
%  (4) both differential operators are integers
%  (5) there is no weight function
%  (6) RNG is equal to the range of the two bases.

%  Else proceed with the use of the Romberg integration.
  
if  isa_fd(fdobj1)               && isa_fd(fdobj2)          && ...
    strcmp(type1,'bspline')      && strcmp(type2,'bspline') && ...
    basisobj1 == basisobj2       && ...
    isinteger(Lfdobj1)           && isinteger(Lfdobj2)      && ...
    nargin < 6                   && ...
    all(rng == range1)
    deriv1 = getnderiv(Lfdobj1);
    deriv2 = getnderiv(Lfdobj2);
    inprodmat = inprod_bspline(fdobj1, fdobj2, deriv1, deriv2);
    iter = 0;
    return;
end

%  ------------------------------------------------------------
%  Now determine the number of subintervals within which the
%  numerical integration takes.  This is important if either
%  basis is a B-spline basis and has multiple knots at a 
%  break point.
%  ------------------------------------------------------------

% The default case, no multiplicities.  

rngvec = rng;  

%  check for any knot multiplicities in either argument

knotmult = [];

%  check first functional object for knot multiplicities

if isa_fd(fdobj1)
    %  univariate case
    if strcmp(type1,'bspline') 
        % Look for knot multiplicities in first basis
        params1  = getbasispar(basisobj1);
        nparams1 = length(params1);
        for i=2:nparams1
            if params1(i) == params1(i-1)
                knotmult = [knotmult, params1(i)];
            end
        end
    end
elseif isa_bifd(fdobj1)
    %  bivariate case
    if strcmp(stype1,'bspline') 
        % Look for knot multiplicities in first basis
        sparams1   = getbasispar(sbasisobj1);
        nsparams1  = length(sparams1);
        for i=2:nsparams1
            if sparams1(i) == sparams1(i-1)
                knotmult = [knotmult, sparams1(i)];
            end
        end
    end
    if strcmp(ttype1,'bspline') 
        % Look for knot multiplicities in first basis
        tparams1   = getbasispar(tbasisobj1);
        ntparams1  = length(tparams1);
        for i=2:ntparams1
            if tparams1(i) == tparams1(i-1)
                knotmult = [knotmult, tparams1(i)];
            end
        end
    end
else
    %  cell case
    basisobj11 = getbasis(fdobj1{1});
    type11 = getbasistype(basisobj11);
    if strcmp(type11,'bspline') 
        % Look for knot multiplicities in first basis
        params11  = getbasispar(basisobj11);
        nparams11 = length(params11);
        for i=2:nparams11
            if params11(i) == params11(i-1)
                knotmult = [knotmult, params11(i)];
            end
        end
    end
    basisobj12 = getbasis(fdobj1{2});
    type12 = getbasistype(basisobj12);
    if strcmp(type12,'bspline') 
        % Look for knot multiplicities in first basis
        params12  = getbasispar(basisobj12);
        nparams12 = length(params12);
        for i=2:nparams12
            if params12(i) == params12(i-1)
                knotmult = [knotmult, params12(i)];
            end
        end
    end
end

%  check second functional object for knot multiplicities

if isa_fd(fdobj2)
    %  univariate case
    if strcmp(type2,'bspline') 
        % Look for knot multiplicities in first basis
        params2  = getbasispar(basisobj2);
        nparams2 = length(params2);
        for i=2:nparams2
            if params2(i) == params2(i-1)
                knotmult = [knotmult, params2(i)];
            end
        end
    end
elseif isa_bifd(fdobj2)
    %  bivariate case
    if strcmp(stype2,'bspline') 
        % Look for knot multiplicities in first basis
        sparams2   = getbasispar(sbasisobj2);
        nsparams2  = length(sparams2);
        for i=2:nsparams2
            if sparams2(i) == sparams2(i-1)
                knotmult = [knotmult, sparams2(i)];
            end
        end
    end
    if strcmp(ttype2,'bspline') 
        % Look for knot multiplicities in first basis
        tparams2   = getbasispar(tbasisobj2);
        ntparams2  = length(tparams2);
        for i=2:ntparams2
            if tparams2(i) == tparams2(i-1)
                knotmult = [knotmult, tparams2(i)];
            end
        end
    end
else
    %  cell case
    basisobj21 = getbasis(fdobj2{1});
    type21 = getbasistype(basisobj21);
    if strcmp(type21,'bspline') 
        % Look for knot multiplicities in first basis
        params21  = getbasispar(basisobj21);
        nparams21 = length(params21);
        for i=2:nparams21
            if params21(i) == params21(i-1)
                knotmult = [knotmult, params21(i)];
            end
        end
    end
    basisobj22 = getbasis(fdobj2{2});
    type22 = getbasistype(basisobj22);
    if strcmp(type22,'bspline') 
        % Look for knot multiplicities in first basis
        params22  = getbasispar(basisobj22);
        nparams22 = length(params22);
        for i=2:nparams22
            if params22(i) == params22(i-1)
                knotmult = [knotmult, params22(i)];
            end
        end
    end
end

%  Modify RNGVEC defining subinvervals if there are any
%  knot multiplicities.

if length(knotmult) > 0
    knotmult = sort(unique(knotmult));
    knotmult = knotmult(knotmult > rng(1) & knotmult < rng(2));
    rngvec = [rng(1), knotmult, rng(2)];
end

inprodmat = zeros(nrep1, nrep2);

%  -----------------------------------------------------------------
%                   loop through sub-intervals
%  -----------------------------------------------------------------

nrng = length(rngvec);
for irng = 2:nrng
    rngi = [rngvec(irng-1),rngvec(irng)];
    %  change range so as to avoid being exactly on
    %  multiple knot values
    if irng > 2
        rngi(1) = rngi(1) + 1e-10;
    end
    if irng < nrng
        rngi(2) = rngi(2) - 1e-10;
    end
    
    %  set up first iteration
    
    iter  = 1;
    width = rngi(2) - rngi(1);
    JMAXP = JMAX + 1;
    h     = ones(JMAXP,1);
    h(2)  = 0.25;
    s = reshape(zeros(JMAXP*nrep1*nrep2,1),[JMAXP,nrep1,nrep2]);
    %  The first iteration uses just the endpoints.  Evaluate
    %  the objects at these endpoints.
    if isa_fd(fdobj1)
        fx1 = eval_fd(rngi, fdobj1, Lfdobj1);
    elseif isa_bifd(fdobj1)
        fx1 = evaldiag_bifd(rngi, fdobj1, 0, Lfdobj1);
    else
        fx11 = eval_fd(rngi, fdobj1{1}, 0);
        fx12 = eval_fd(rngi, fdobj1{2}, Lfdobj1);
        if     coefd11(2) == 1 && coefd12(2) > 1
            fx1 = (fx11*ones(1,coefd12(2))).*fx12;
        elseif coefd11(2) > 1 && coefd12(2) == 1
            fx1 = fx11.*(fx12*ones(1,coefd11(2)));
        else
            fx1  = fx11.*fx12;
        end
    end
    if isa_fd(fdobj2)
        fx2 = eval_fd(rngi, fdobj2, Lfdobj2);
    elseif isa_bifd(fdobj2)
        fx2 = evaldiag_bifd(rngi, fdobj2, 0, Lfdobj2);
    else
        fx21 = eval_fd(rngi, fdobj2{1}, 0);
        fx22 = eval_fd(rngi, fdobj2{2}, Lfdobj2);
        if     coefd21(2) == 1 && coefd22(2) > 1
            fx2 = (fx21*ones(1,coefd22(2))).*fx22;
        elseif coefd21(2) > 1 && coefd22(2) == 1
            fx2 = fx21.*(fx22*ones(1,coefd21(2)));
        else
            fx2  = fx21.*fx22;
        end
    end
    %  multiply by values of weight function if necessary
    if ~isnumeric(wtfd)
        wtd = eval_fd(wtfd, rngi);
        fx2 = (wtd * ones(1,nrep2)) .* fx2;
    end
    s(1,:,:)  = width.*(fx1' * fx2)./2;
    tnm = 0.5;
    
    %  now iterate to convergence
    
    for iter = 2:JMAX
        tnm = tnm.*2;
        del = width./tnm;
        x   = rngi(1)+del/2:del:rngi(2);
        if isa_fd(fdobj1)
            fx1 = eval_fd(x, fdobj1, Lfdobj1);
        elseif isa_bifd(fdobj1)
            fx1 = evaldiag_bifd(x, fdobj1, 0, Lfdobj1);
        else
            fx11 = eval_fd(x, fdobj1{1}, 0);
            fx12 = eval_fd(x, fdobj1{2}, Lfdobj1);
            if     coefd11(2) == 1 && coefd12(2) > 1
                fx1 = (fx11*ones(1,coefd12(2))).*fx12;
            elseif coefd11(2) > 1 && coefd12(2) == 1
                fx1 = fx11.*(fx12*ones(1,coefd11(2)));
            else
                fx1  = fx11.*fx12;
            end
        end
        if isa_fd(fdobj2)
            fx2 = eval_fd(x, fdobj2, Lfdobj2);
        elseif isa_bifd(fdobj2)
            fx2 = evaldiag_bifd(x, fdobj2, 0, Lfdobj2);
        else
            fx21 = eval_fd(x, fdobj2{1}, 0);
            fx22 = eval_fd(x, fdobj2{2}, Lfdobj2);
            if     coefd21(2) == 1 && coefd22(2) > 1
                fx2 = (fx21*ones(1,coefd22(2))).*fx22;
            elseif coefd21(2) > 1 && coefd22(2) == 1
                fx2 = fx21.*(fx22*ones(1,coefd21(2)));
            else
                fx2  = fx21.*fx22;
            end
        end
        %  multiply by values of weight function if necessary
        if ~isnumeric(wtfd)
            wtd = eval_fd(wtfd, x);
            fx2 = (wtd * ones(1,nrep2)) .* fx2;
        end
        chs = reshape(width.*(fx1' * fx2)./tnm,[1,nrep1,nrep2]);
        s(iter,:,:) = (s(iter-1,:,:) + chs)./2;
        if iter >= 5
            ind = (iter-4):iter;
            ya = s(ind,:,:);
            xa = h(ind);
            absxa = abs(xa);
            [absxamin, ns] = min(absxa);
            cs = ya;
            ds = ya;
            y  = squeeze(ya(ns,:,:));
            ns = ns - 1;
            for m = 1:4
                for i = 1:(5-m)
                    ho      = xa(i);
                    hp      = xa(i+m);
                    w       = (cs(i+1,:,:) - ds(i,:,:))./(ho - hp);
                    ds(i,:,:) = hp.*w;
                    cs(i,:,:) = ho.*w;
                end
                if 2*ns < 5-m
                    dy = squeeze(cs(ns+1,:,:));
                else
                    dy = squeeze(ds(ns,:,:));
                    ns = ns - 1;
                end
                y = y + dy;
            end
            ss     = reshape(y, nrep1, nrep2);
            errval = max(max(abs(dy)));
            ssqval = max(max(abs(ss)));
            if all(ssqval > 0)
                crit = errval./ssqval;
            else
                crit = errval;
            end
            if crit < EPS && iter >= JMIN
                break
            end
        end
        s(iter+1,:,:) = s(iter,:,:);
        h(iter+1)   = 0.25.*h(iter);
        if iter == JMAX
            warning('Wid:converge','Failure to converge.');
        end
    end
    
    inprodmat = inprodmat + ss;
    
end

