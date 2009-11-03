function linmodstr = linmod_new(xfd, yfd, wtvec, xfdPar, yfdPar, zmatrnk)
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

if nargin < 6, zmatrnk = [];          end

if isa(xfd, 'fd') && isa(yfd, 'fd')
   
    if nargin < 5, yfdPar = fdPar(getbasis(yfd));  end
    if nargin < 4, xfdPar = fdPar(getbasis(xfd));  end
    
    xlambda = getlambda(xfdPar);
    ylambda = getlambda(yfdPar);
    xLfdobj = getLfd(xfdPar);
    yLfdobj = getLfd(yfdPar);
    
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
    if nargin < 3 || isempty(wtvec), wtvec = ones(ncurves,1); end
    
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


