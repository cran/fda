function linmodstr = linmod(xfd, yfd, wtvec, ...
                   xLfdobj, yLfdobj, xlambda, ylambda, zmatrnk)
%  LINMOD  Fits a functional linear model. t
%  The model consists of a constant term
%  plus either a conventional independent variable matrix or a single
%  functional independent variable.  The modeling problem may
%  be functional either in terms of the independent variable or in
%  terms of the dependent variable, but at least one variable must
%  be functional.
%  Smoothing is controlled by two parameters XLAMBDA and YLAMBDA,
%  corresponding to the independent and dependent functional
%  variables, respectively.
%  Arguments:
%  XFD     ... a data struct for the independent variable that may
%              be either of 'fd' class or a matrix
%  YFD     ... a data struct for the   dependent variable that may
%              be either of 'fd' class or a matrix
%  WTVEC   ... a vector of weights
%  XLFDOBJ ... either an integer defining a derivative, or a structure
%                defining a linear differential operator to be applied
%                to penalize roughness in the argument for xfd
%  YFDOBJ ... as XLfdobj, but for roughness in the argument for yfd
%  XLAMBDA ... a smoothing parameter for the independent variable arg.
%  YLAMBDA ... a smoothing parameter for the   dependent variable arg.
%  ZMATRNK ... actual rank of independent variable matrix for the
%              functional DV/multivariate IV case
%  Returns:  a struct object LINMODSTR with fields
%  ALPHA ... intercept term (can be either vector or fd)
%  REG   ... a functional data struct for the regression function
%  YHAT  ... the approximation to y, which can be a vector or fd

%  Last modified 2 December 2006

%  set up some default argument values

if nargin < 8, zmatrnk = [];          end
if nargin < 7, ylambda = 0;           end
if nargin < 6, xlambda = 0;           end
if nargin < 5, xLfdobj = int2Lfd(2);  end
if nargin < 4, yLfdobj = int2Lfd(2);  end

%  check XLFDOBJ and YLFDOBJ

xLfdobj = int2Lfd(xLfdobj);
yLfdobj = int2Lfd(yLfdobj);

%  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%              The multivariate IV and functional DV case
%  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if isa(yfd, 'fd') && isa(xfd, 'double')
    
    if ~isa_fd(yfd)
        error ('Argument yfd is not a functional data object.');
    end
    
    if ndims(xfd) == 2
        zmat  = xfd;
        fdobj = yfd;
    else
        error (['First argument not a two-dimensional array', ...
                'when second argument is functional']);
    end
    
    coef   = getcoef(fdobj);
    coefd  = size(coef);
    ndim   = length(coefd);
    if (ndim < 2)
        error('Linear modeling impossible with 1 replication');
    end
    
    [ncurves,p] = size(zmat);
    if nargin < 3, wtvec = ones(ncurves,1); end
    
    if(ndim == 3)
        nvar = coefd(3);
    else
        nvar = 1;
    end
    
    basisobj = getbasis(fdobj);
    nbasis   = getnbasis(basisobj);
    
    if nargin < 8
        zmatrnk = p;
    end
    
    %  check weight vector
    
    if nargin < 3
        wtvec = ones(ncurves,1);
    end
    
    if length(wtvec) ~= ncurves
        error('WTVEC of wrong length');
    end
    
    rangewt = [min(wtvec), max(wtvec)];
    
    if rangewt(1) < 0
        error('WTVEC must not contain negative values.');
    end
    
    if min(wtvec) <= 0
        error('All values of WTVEC must be positive.');
    end
    
    if coefd(2) ~= ncurves
        error('Number of rows of ZMAT must equal number of replications');
    end
    
    zd = size(zmat);
    p  = zd(2);
    
    if zmatrnk > p
        error('Specified rank of ZMAT must not be greater than no. columns.');
    end
    
    if nvar > 1
        bcoef = zeros(nbasis,p,nvar);
    else
        bcoef = zeros(nbasis,p);
    end
    
    %  rescale if weights not constant
    
    if rangewt(1) ~= rangewt(2)
        rootwt = sqrt(wtvec);
        if nvar == 1
            zmat = zmat .* (rootwt * ones(1,p));
            coef = coef .* (ones(nbasis,1) * rootwt');
        else
            for j=1:nvar
                zmat(:,:,j) = zmat(:,:,j)  .* (rootwt * ones(1,p));
                coef(:,:,j) = coef(:,:,j)' .* (ones(nbasis,1) * rootwt');
            end
        end
    end
    
    if ylambda < 0
        ylambda = 0;
    end
    
    nyderiv = getnderiv(yLfdobj);
    if ylambda > 0 && nyderiv < 0
        error('Order of derivative must be nonnegative.');
    end
    
    if zmatrnk < p
        [zmatu,zmatd,zmatv] = svd(zmat,0);
        zmatd = diag(zmatd);
        if zmatd(zmatrnk) <= 0
            error('ZMAT is not of specified column rank');
        end
        index  = 1:zmatrnk;
        zginvt = zmatu(:,index) * diag(1/zmatd(index)) * zmatv(:,index)';
        if nvar == 1
            bcoef = coef * zginvt;
        else
            for j = 1:nvar
                bcoef(:,:,j) = coef(:,:,j) * zginvt;
            end
        end
    else
        if nvar == 1
            bcoef = (zmat\coef')';
        else
            for j = 1:nvar
                bcoef(:,:,j) = (zmat\coef(:,:,j)')';
            end
        end
    end
    
    yhatcoef = bcoef * zmat';
    
    fdnames    = getnames(yfd);
    fdnames{2} = 'Reg. Coefficients';
    regfd  = fd(bcoef,    basisobj, fdnames);
    fdnames    = getnames(yfd);
    fdnames{2} = ['Predicted ', fdnames{2}];
    yhatfd = fd(yhatcoef, basisobj, fdnames);
    
    linmodstr.alpha = 0;
    linmodstr.reg   = regfd;
    linmodstr.yhat  = yhatfd;
    
    return;
    
end

%  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             The functional IV and multivariate DV case
%  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if isa(xfd, 'fd') && isa(yfd, 'double')
    
    if ~isa_fd(xfd)
        error ('Argument xfd is not a functional data object.');
    end
    
    sizey = size(yfd);
    if length(sizey) == 2
        ymat  = yfd;
        fdobj = xfd;
    else
        error (['Second argument not a functional data object', ...
                'when first argument is functional.']);
    end
    
    coef  = getcoef(fdobj);
    coefd = size(coef);
    ndim  = length(coefd);
    if ndim < 2
        error('Linear modeling impossible with 1 replication');
    end
    if ndim == 3
        error('This version cannot accommodate multiple functional IVs');
    end
    ncurves = coefd(2);
    if nargin < 3, wtvec = ones(ncurves,1); end
    
    %  check weight vector
    
    if length(wtvec) ~= ncurves
        error('WTVEC of wrong length');
    end
    rangewt = [min(wtvec), max(wtvec)];
    if (rangewt(1) <= 0)
        error('WTVEC must not contain negative values.');
    end
    if (min(wtvec) <= 0)
        error('All values of WTVEC must be positive.');
    end
    
    basisobj = getbasis(fdobj);
    nbasis   = getnbasis(basisobj);
    if sizey(1) ~= ncurves
        error('Number of rows of YMAT must equal number of replications');
    end
    
    if nargin < 6
        xlambda = 0;
    end
    
    if nargin < 4
        xLfdobj = int2Lfd(2);
    end
    
    if xlambda < 0
        error ('Value of XLAMBDA was negative.');
    end
    
    jmat = inprod(basisobj, basisobj);
    zmat = [ones(1,ncurves); jmat * coef]';
    
    index = 2:(nbasis+1);
    
    if xlambda <= 0
        %  no smoothing required, do ordinary least squares
        temp = zmat\ymat;
    else
        %  smoothing required
        kmat = zeros(nbasis+1,nbasis+1);
        kmat(index,index) = inprod(basisobj, basisobj, xLfdobj, xLfdobj);
        Cmat = zmat' * zmat + xlambda .* kmat;
        Dmat = zmat' * ymat;
        temp = symsolve(Cmat,Dmat);
    end
    bcoef = temp(index,:);
    alpha = temp(1,:);
    yhat  = zmat * temp;
    
    regfdnames = getnames(xfd);
    regfdnames{3} = 'Reg. Coefficient';
    regfd = fd(bcoef, basisobj, regfdnames);
    
    linmodstr.alpha = alpha;
    linmodstr.reg   = regfd;
    linmodstr.yhat  = yhat;
    
    return;
    
end

%  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             The functional IV and functional DV case
%  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if isa(xfd, 'fd') && isa(yfd, 'fd')
    
    if ~isa_fd(xfd)
        error ('Argument xfd is not a functional data object.');
    end
    
    if ~isa_fd(yfd)
        error ('Argument yfd is not a functional data object.');
    end
    
    coefx   = getcoef(xfd);
    coefy   = getcoef(yfd);
    coefdx  = size(coefx);
    coefdy  = size(coefy);
    ndimx   = length(coefdx);
    ndimy   = length(coefdy);
    
    if ndimx < 2
        error('Linear modeling impossible with 1 replication');
    end
    
    if ndimx == 3
        error('This version cannot accommodate multiple functional IVs');
    end
    
    ncurves = coefdx(2);
    if coefdy(2) ~= ncurves
        error ('Numbers of observations in first two arguments do not match.');
    end
    onen = ones(1,ncurves);
    if nargin < 3, wtvec = ones(ncurves,1); end
    
    rangewt = [min(wtvec),max(wtvec)];
    if rangewt(1) < 0
        error('WTVEC must not contain negative values.');
    end
    if rangewt(1) ~= rangewt(2)
        error('Unequal weight option not enabled');
    end
    
    xbasisobj = getbasis(xfd);
    ybasisobj = getbasis(yfd);
    
    nbasisx = coefdx(1);
    nbasisy = coefdy(1);
    
    if length(wtvec) ~= ncurves
        error('WTVEC of wrong length');
    end
    if min(wtvec) <= 0
        error('All values of WTVEC must be positive.');
    end
    
    if nargin < 7
        ylambda = 0;
    end
    
    if nargin < 6
        xlambda = 0;
    end
    
    if nargin < 5
        yLfdobj = int2Lfd(2);
    end
    
    if nargin < 4
        xLfdobj = int2Lfd(2);
    end
    
    if xlambda < 0
        warning ('Wid:negative', ...
            'Value of LAMBDA was negative, and 0 used instead.');
    end
    
    jmatx   = inprod(xbasisobj, xbasisobj);
    penmatx = inprod(xbasisobj, xbasisobj, xLfdobj, xLfdobj);
    if ndimx == 2
        zmatx   = [onen; jmatx * coefx]';
    else
        zmatx   = [onen; jmatx * coefx(:,:,1)]';
    end
    
    jmaty   = inprod(ybasisobj, ybasisobj);
    penmaty = inprod(ybasisobj, ybasisobj, yLfdobj, yLfdobj);
    
    index = 2:(nbasisx+1);
    kmatx = zeros(nbasisx+1,nbasisx+1);
    kmatx(index,index) = penmatx;
    
    tempx = inv(zmatx' * zmatx + xlambda*kmatx);
    tempy = inv(jmaty          + ylambda*penmaty);
    if ndimy == 2
        gmat = tempx * zmatx' * coefy'        * jmaty * tempy;
    else
        gmat = tempx * zmatx' * coefy(:,:,1)' * jmaty * tempy;
    end
    yhatcoef  = (zmatx * gmat)';
    
    bcoef = zeros(nbasisx,nbasisy,1,1);
    bcoef(:,:,1,1) = gmat(index,:);
    acoef = gmat(1,:)';
    
    %  functional structure for the alpha function
    
    alphafdnames = getnames(yfd);
    alphafdnames{3} = 'Intercept';
    alphafd = fd(acoef, ybasisobj, alphafdnames);
    
    %  bi-functional structure for the beta function
    
    regfdnames = getnames(xfd);
    regfdnames{3} = 'Reg. Coefficient';
    regfd = bifd(bcoef, xbasisobj, ybasisobj, regfdnames);
    
    %  functional data structure for the yhat functions
    
    yhatfd = fd(yhatcoef, ybasisobj, getnames(yfd));
    
    linmodstr.alpha = alphafd;
    linmodstr.reg   = regfd;
    linmodstr.yhat  = yhatfd;
    
    return;
end


