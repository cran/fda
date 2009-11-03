function [bwtcell, awtcell, resfdcell] = ...
    pda_fd(xfdcell, bwtcell, awtcell, ufdcell, norder, nfine)

%PDA_FD computes the basis function expansions of the
%  estimates of the coefficient functions a_k(t) and b_j(t) 
%  in the possibly nonhomogeneous linear differential operator
%
%    Lx(t) = b_0(t)x(t) + b_1(t)Dx(t) + ... + b_{m-1}D^{m-1}x(t) + D^m x(t)
%          - a_1(t)u_1(t) - ... - a_k(t)u_K(t) 
%  of order m = NORDER that minimizes in a least squares sense the 
%  residual functions f(t) = Lx(t).  
%
%  If NORDER = 0, PDACELL fits the varying coefficient or pointwise
%  linear model using the functions x(t) as dependent variables and
%  the forcing functions u(t) as independent variables.  In this case,
%  there must be at least one forcing function.
%
%  The functions x(t) are in functional data object XFDOBJ.
%  The forcing functions u_k(t) are in functional data object UFDOBJ.
%  The coefficient functions for u_k(t) and x(t) are expanded in terms 
%  of basis functions specified in AWTCELL and BWTCELL, respectively.
%
%  The functions u_k(t) and x(t) are assumed to be vector valued 
%    of dimension J. 
%  That is, the differential equation can be a system of J equations 
%    rather than a single equation.   
%  Each coefficient function b_j(t) is matrix valued, with a column
%    for each scalar function in the system.  

%  Arguments:
%  XFDCELL   ...  cell array of functional data objects for the 
%                 functions whose derivatives define the DIFE
%                 dimensions are J and 1
%  BWTCELL   ...  cell array of weight function specs for functions x
%                 dimensions are J, J and NORDER. 
%                 Each cell contains an fdPar object.
%  AWTCELL   ...  cell array of weight function specs for u-variables
%                 dimensions are J and K.   
%                 Each cell contains an fdPar object.
%  UFDCELL   ...  cell array of independent variables or u-variables
%                 dimensions are J and K
%  NORDER    ...  order of the linear differential operator, that is, 
%                 the order of the highest derivative.
%  NFINE     ...  number of sampling points for numerical integration
%                 Warning:  NFINE can make quite a difference.  If all 
%                 functions are relatively smooth, a value in the low
%                 hundreds may be fine, but if there are sharp break
%                 points, as there would be in very low error situations
%                 with abrupt changes in inputs, a thousand for more
%                 sampling points may be required.  Ultimately adaptive
%                 quadrature is needed here, as elsewhere in FDA.    

%  Returns:
%  BFDCELL   ...  cell array of weights for x functions
%                 dimensions are J, J and NORDER
%  AFDCELL   ...  cell array of weights for u-variables
%                 dimension J and K
%  RESFDCELL ...  FD object for residual functions.

%  last modified 2 November 2008

if nargin < 2
    error('There are less than three arguments.');
end

%   set some default arguments

if nargin < 6, nfine   = 501;  end
if nargin < 4, ufdcell = {};   end
if nargin < 3, awtcell = {};   end

%  get and check NORDER

if nargin < 5, 
    bwtdims = size(bwtcell);
    if length(bwtdims) > 3
        error('BWTCELL has more than three dimensions.');
    end
    if all(bwtdims) == 1 || length(bwtdims) == 2
        norder = 1;
    else
        norder = bwtdims(3);
    end
end

norder = floor(norder);
if norder < 0, error('NORDER is negative.');  end
nordp1 = norder + 1;

%  check dimensions of cells

nvar = size(xfdcell,1);
if size(xfdcell,2) ~= 1
    error('XFDCELL has more than one column.');
end

%  ----------------------------------------------------------------
%     For efficiency, there are two versions of this code:
%     one for a single variable, and another for multiple variables.
%  ----------------------------------------------------------------

%  ----------------------------------------------------------------
%                   Single variable case
%  ----------------------------------------------------------------

if nvar == 1
    
    %  check the dimensions of UFDCELL and AWTCELL
    
    if isempty(ufdcell) || isempty(awtcell)
        nforce  = 0;
        afdcell = {};
    else
        nforce = length(ufdcell);
        urange = getbasisrange(getbasis(ufdcell{1}));
        if length(awtcell) ~= nforce
            error('The length of AWTCELL is incorrect.');
        end
        for k=1:nforce
            if ~isa_fd(ufdcell{k})
                error(['UFDCELL{',num2str(k), ...
                       '} is not a functional data object.'])
            end
            if ~isa_fdPar(awtcell{k})
                error(['AWTCELL{',num2str(k), ...
                       '} is not a functional parameter object.'])
            end
        end
    end
    
    %  check to see if there is anything to estimate
    
    if norder == 0 && nforce == 0
        error('There are no coefficient functions to estimate.');
    end
    
    %  check the dimensions of BWTCELL
    
    if norder == 0
        bfdcell = {};
    end
    if length(bwtcell) ~= norder
        error('The dimensions of BWTCELL are incorrect.');
    end
    
    %  check XFDCELL and extract NCURVE and XRANGE
    
    if ~isa_fd(xfdcell{1})
        error('XFDCELL{1} is not a functional data object.');
    end 
    xbasis     = getbasis(xfdcell{1});
    xrange     = getbasisrange(xbasis);
    nxbasis    = getnbasis(xbasis);
    ncurve     = size(getcoef(xfdcell{1}),2);
    bfdnames   = getnames(xfdcell{1});
    resfdnames = bfdnames;
    
    nbasmax = nxbasis;  %  This will be the maximum number of basis functions
    
    %  check UFDCELL and extract URANGE
    %  Note:  XRANGE and URANGE are required to be identical 
    %         in this version
    
    if nforce > 0
        for k=1:nforce
            if k == 1
                urange   = getbasisrange(getbasis(ufdcell{k}));
                %  check that URANGE is equal to XRANGE
                if any(urange ~= xrange)
                    error('XRANGE and URANGE are not identical');
                end
            else
                if any(getbasisrange(getbasis(ufdcell{k})) ~= urange)
                    error('Ranges are incompatible for UFDCELL.');
                end
            end
        end
        
        %  check AWTCELL and extract the max. no. basis fns.
        
        for k=1:nforce
            afdi = getfd(awtcell{k});
            if ~isa_fd(afdi)
                error('AFDI is not a functional data object.');
            end  
            basisi = getbasis(afdi);
            if any(getbasisrange(basisi) ~= urange)
                error('Ranges are incompatible for AWTCELL.');
            end
            nbasi   = getnbasis(basisi);
            nbasmax = max([nbasmax,nbasi]);
        end
        
    end
    
    %  check BWTCELL
    
    if norder > 0
        for j=1:norder
            bfdPar12 = bwtcell{j};
            if ~isa_fdPar(bfdPar12)
              error('BWTCELL{1} is not a functional parameter object.');
            end
            bfd12  = getfd(bfdPar12);
            basisi = getbasis(bfd12);
            if any(getbasisrange(basisi) ~= xrange)
                error('Ranges are incompatible in BWTCELL.');
            end
            nbasi   = getnbasis(basisi);
            nbasmax = max([nbasmax,nbasi]);
        end  
    end
    
    %  Set up sampling values to be used in numerical integration
    %    and set up matrix of basis values.  The number of sampling
    %  NFINE is here set to a usually workable value if too small.
    
    if nargin < 6, nfine = 501;  end
    if nfine  < 5*nbasmax, nfine = 5*nbasmax;  end
    
    deltax = (xrange(2)-xrange(1))/(nfine-1);
    tx     = (xrange(1):deltax:xrange(2))';
    deltau = deltax;
    tu     = tx;
    
    %  set up  YARRAY to hold values of x functions and 
    %   their derivatives
    
    yarray = zeros([nfine,ncurve,nordp1]);
    for j=1:nordp1
        yarray(:,:,j) = eval_fd(tx, xfdcell{1}, j-1);
    end
    
    %  set up  UARRAY to hold values of u functions
    
    if nforce > 0
        uarray = zeros([nfine,ncurve,nforce]);
        for k=1:nforce
            temp = eval_fd(tu, ufdcell{k});
            uarray(:,:,k) = temp;
        end
    end
    
    %  set up array YPROD to hold mean of products of values in YARRAY
    
    yprod = zeros([nfine,nordp1,nordp1]);
    for j1=1:nordp1
        for j2=1:j1;
            if ncurve == 1
                yprodval = squeeze(yarray(:,1,j1)).* ...
                           squeeze(yarray(:,1,j2));
            else
                yprodval = mean(squeeze(yarray(:,:,j1)).* ...
                    squeeze(yarray(:,:,j2)),2);
            end
            yprod(:,j1,j2) = yprodval;
            yprod(:,j2,j1) = yprodval;
        end
    end
    
    %  set up array YUPROD to hold mean of u-variables u times 
    %    x functions and their derivatives
    
    onesncurve = ones(1,ncurve);
    if nforce > 0
        yuprod = zeros(nfine, nforce, nordp1);
        for k=1:nforce
            for j1=1:nordp1
                if ncurve == 1
                    yuprodval = yarray(:,1,j1).*uarray(:,k);
                else
                    yuprodval = ...
                        mean(squeeze(yarray(:,:,j1).* ...
                        uarray(:,:,k)),2);
                end
                yuprod(:,k,j1) = yuprodval;
            end
        end
    end
    
    clear yarray
    
    %  set up array UPROD to hold mean of products of u-variables u 
    
    if nforce > 0
        uprod = zeros(nfine, nforce, nforce);
        for k=1:nforce
            for ju=1:k
                if ncurve == 1
                    uprodval = uarray(:,1,k).*uarray(:,1,ju);
                else
                    uprodval = mean(squeeze(uarray(:,:,k).* ...
                        uarray(:,:,ju)),2);           
                end
                uprod(:,k,ju) = uprodval;
                uprod(:,ju,k) = uprodval;
            end
        end
        clear uarray
    end
    
    %  set up an index array and some arrays of 1's
    
    onesn = ones(nfine,1);
    
    %  set up array to hold coefficients for basis expansions
    
    if nforce > 0
        aarray = zeros(nfine,nforce);  
    else
        aarray = [];
    end
    
    if norder > 0
        barray = zeros(nfine,norder);  
    else
        barray = [];
    end
    
    
    %  -------  beginning of loop through variables  --------------
    
    %  get number of coefficients to be estimated for this equation
    
    % loop through u-variables
    neqns  = 0;
    for k = 1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark) 
            neqns = neqns + getnbasis(getbasis(getfd(afdPark)));
        end
    end
    % loop through x functions and their derivatives
    for j2=1:norder
        bfdParij = bwtcell{j2};
        if getestimate(bfdParij)
            neqns = neqns + getnbasis(getbasis(getfd(bfdParij)));
        end
    end
    if neqns < 1
        error('Number of equations to solve is not positive.');
    end
    
    %  set up coefficient array and right side array for 
    %     linear equation
    
    cmat   = zeros(neqns, neqns);
    dmat   = zeros(neqns, 1);
    
    %  evaluate default weight functions for this variable
    
    for k=1:nforce
        afdPark = awtcell{k};
        aarray(:,k) = eval_fd(tu, getfd(afdPark));
    end
    
    for j=1:norder
        bfdParij = bwtcell{j};
        barray(:,j) = squeeze(eval_fd(tx, getfd(bfdParij)));
    end
    
    %  loop through equations, 
    %    corresponding to rows for CMAT and DMAT
    
    %  loop through equations for u-variables
    
    mi12 = 0;
    for k1 = 1:nforce
        afdPark1   = awtcell{k1};
        if getestimate(afdPark1)
            abasisk1    = getbasis(getfd(afdPark1));
            abasismatk1 = getbasismatrix(tu, abasisk1);
            mi11 = mi12 + 1;
            mi12 = mi12 + getnbasis(abasisk1);
            indexi1 = mi11:mi12;
            %  DMAT entry for u-variable
            weightk1 = -yuprod(:,k1,nordp1);
            dmat(indexi1) = trapzmat(abasismatk1, onesn, deltax, ...
                                     weightk1);
            %  add terms corresponding to x-derivative weights
            %  that are not estimated
            for j=1:norder;
                bfdParij   = bwtcell{j};
                if ~getestimate(bfdParij)
                    weightij = -yuprod(:,k1,j);
                    dmat(indexi1) = dmat(indexi1) + ...
                        trapzmat(abasismatk1, barray(:,j), deltax, ...
                                 weightij);
                end
            end
            %  loop through weight functions to be estimated,
            %    corresponding to columns for CMAT
            %  begin with u-variable -- u-variable pairs
            mi22 = 0;
            for k2=1:nforce
                afdPark2   = awtcell{k2};
                %  if estimation required, modify coefficient matrix
                if getestimate(afdPark2)
                    abasisk2    = getbasis(getfd(afdPark2));
                    abasismatk2 = getbasismatrix(tu, abasisk2);
                    weightk2    = uprod(:,k1,k2);
                    Cprod       = ...
                        trapzmat(abasismatk1, abasismatk2, deltau, ...
                                 weightk2);
                    mi21 = mi22 + 1;
                    mi22 = mi22 + getnbasis(abasisk2);
                    indexi2 = mi21:mi22;
                    cmat(indexi1,indexi2) = Cprod;
                end
            end
            %  remaining columns: 
            %    loop through u-variable -- x-derivative pairs
            mij22 = mi22;
            for j2=1:norder;
                bfdParij2   = bwtcell{j2};
                %  if estimation required, modify coefficient matrix
                if getestimate(bfdParij2)
                    bbasisij2    = getbasis(getfd(bfdParij2));
                    bbasismatij2 = getbasismatrix(tx, bbasisij2);
                    weightij12   = -yuprod(:,k1,j2);
                    Cprod = ...
                        trapzmat(abasismatk1, bbasismatij2, deltax, ...
                                 weightij12);
                    mij21 = mij22 + 1;
                    mij22 = mij22 + getnbasis(bbasisij2);
                    indexij2  = mij21:mij22;
                    cmat(indexi1,indexij2) = Cprod;
                end
            end
            %  add roughness penalty matrix to diagonal entries
            lambdakl = getlambda(afdPark1);
            if lambdakl > 0.0
                Lfdobj = getLfd(afdPark1);
                penmat = lambdakl.*eval_penalty(abasisk1, Lfdobj);
                cmat(indexi1,indexi1) = cmat(indexi1,indexi1) + penmat;
            end
        end
    end
    
    %  loop through equations for x-derivatives
    
    mij12 = mi12;
    for j1=1:norder
        bfdParij1 = bwtcell{j1};
        if getestimate(bfdParij1)
            bbasisij1    = getbasis(getfd(bfdParij1));
            bbasismatij1 = getbasismatrix(tx, bbasisij1);
            mij11 = mij12 + 1;
            mij12 = mij12 + getnbasis(bbasisij1);
            indexij1 = mij11:mij12;
            %  DMAT entry 
            weightij1 = yprod(:,j1,nordp1);
            dmat(indexij1) = ...
                trapzmat(bbasismatij1,onesn,deltax,weightij1);
            %  add terms corresponding to forcing functions with
            %  unestimated coefficients
            for k=1:nforce
                afdPark   = awtcell{k};
                if ~getestimate(afdPark)
                    weightijk = yuprod(:,k,j1);
                    dmat(indexij1) = dmat(indexij1) + ...
                        trapzmat(bbasismatij1, aarray(:,k), deltax, ...
                                 weightijk);
                end
            end
            %  first columns of CMAT: x-derivative -- u-variable entries
            mi22 = 0;
            for k2=1:nforce
                afdPark2   = awtcell{k2};
                %  if estimation required, modify coefficient matrix
                if getestimate(afdPark2)
                    abasisk2    = getbasis(getfd(afdPark2));
                    abasismatk2 = getbasismatrix(tx, abasisk2);
                    weightk2    = -yuprod(:,k2,j1);
                    Cprod = ...
                        trapzmat(bbasismatij1, abasismatk2, deltax, ...
                                 weightk2);
                    mi21 = mi22 + 1;
                    mi22 = mi22 + getnbasis(abasisk2);
                    indexi2 = mi21:mi22;
                    cmat(indexij1,indexi2) = Cprod;
                end
            end
            %  remaining columns: x-derivative -- x-derivative pairs
            mij22 = mi22;
            for j2=1:norder;
                bfdParij2 = bwtcell{j2};
                bbasisij2    = getbasis(getfd(bfdParij2));
                %  if estimation required, modify coefficient matrix
                if getestimate(bfdParij2)
                    bbasismatij2 = getbasismatrix(tx, bbasisij2);
                    weightij12   = yprod(:,j1,j2);
                    Cprod = ...
                        trapzmat(bbasismatij1, bbasismatij2, deltax, ...
                                 weightij12);
                    mij21 = mij22 + 1;
                    mij22 = mij22 + getnbasis(bbasisij2);
                    indexij2 = mij21:mij22;
                    cmat(indexij1,indexij2) = Cprod;
                end
            end
            %  add roughness penalty matrix to diagonal entries
            lambdaij1 = getlambda(bfdParij1);
            if lambdaij1 > 0.0
                Lfdobj = getLfd(bfdParij1);
                penmat = lambdaij1.*eval_penalty(bbasisij1, Lfdobj);
                cmat(indexij1,indexij1) = cmat(indexij1,indexij1) + ...
                                          penmat;
            end
        end
    end
    
    %  solve for coefficients of basis expansions
    
    dvec = -symsolve(cmat,dmat);
    
    %  set up u-function weight functions
    
    mi2 = 0;
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark)
            afdk = getfd(afdPark);
            mi1 = mi2 + 1;
            mi2 = mi2 + getnbasis(getbasis(afdk));
            indexi = mi1:mi2;
            awtcell{k} = ...
                putfd(afdPark, putcoef(afdk, dvec(indexi)));
        else
            awtcell{k} = afdPark;
        end
    end
    
    %  set up X-function derivative weight functions
    
    mij2 = mi2;
    for j1=1:norder;
        bfdParij = bwtcell{j1};
        if getestimate(bfdParij)
            bfdij = getfd(bfdParij);
            mij1 = mij2 + 1;
            mij2 = mij2 + getnbasis(getbasis(bfdij));
            indexij = mij1:mij2;
            bwtcell{j1} = ...
                putfd(bfdParij, putcoef(bfdij, dvec(indexij)));
        else
            bwtcell{j1} = bfdParij;
        end
    end
    
    %  set up residual cell RESFDCELL
    
    resfdcell     = cell(nvar);
    resfdnames{2} = 'Residual function';
    resfdnames{3} = 'Residual function value';
    
    resbasis = getbasis(xfdcell{1});
    %  initialize with highest order derivative for this variable
    resmat  = eval_fd(tx, xfdcell{1}, norder);
    %  add contributions from weighted u-functions
    for k=1:nforce
        amati    = eval_fd(tu, getfd(awtcell{k}));
        umati    = eval_fd(tu, ufdcell{k});
        resmat   = resmat - (amati*onesncurve).*umati;
    end
    %  add contributions from weighted x-function derivatives
    for j1=1:norder;
        bmatij = eval_fd(tx, getfd(bwtcell{j1}))*onesncurve;
        xmatij = eval_fd(tx, xfdcell{1}, j1-1);
        resmat = resmat + bmatij.*xmatij;
    end
    %  set up the functional data object
    resfdi = smooth_basis(tx, resmat, resbasis);
    resfdi = putnames(resfdi, resfdnames);
    resfdcell{1} = resfdi;

else
    
    %  ----------------------------------------------------------------
    %                   Multiple variable case
    %  ----------------------------------------------------------------
    
    %  check the dimensions of UFDCELL and AWTCELL if there are
    %  
    if isempty(ufdcell) || isempty(awtcell)
        nforce = 0;
        afdcell = {};
    else
        nforce = size(ufdcell,2);
        urange = getbasisrange(getbasis(ufdcell{1}));
        if size(ufdcell,1) ~= nvar
            error(['The number of rows of UFDCELL]', ...
                    ' does not match that of XFDCELL.']);
        end
        if any(size(awtcell) ~= [nvar, nforce])
            error('The dimensions of AWTCELL are incorrect.');
        end
    end
    
    %  check to see if there is anything to estimate
    
    if norder == 0 && nforce == 0
        error('There are no coefficient functions to estimate.');
    end
    
    %  check the dimensions of BWTCELL
    
    if norder == 0
        bfdcell = {};
    end
    if norder == 1
        if any(size(bwtcell) ~= [nvar, nvar])
            error('The dimensions of BWTCELL are incorrect.');
        end
    end
    if norder > 1
        if any(size(bwtcell) ~= [nvar, nvar, norder])
            error('The dimensions of BWTCELL are incorrect.');
        end
    end
    
    %  check XFDCELL and extract NCURVE and XRANGE
    
    for ivar1=1:nvar
        if ~isa_fd(xfdcell{ivar1})
            error(['XFDCELL{',num2str(ivar1), ...
                    '} is not a functional data object.']);
        end 
        if ivar1 == 1
            xrange     = getbasisrange(getbasis(xfdcell{ivar1}));
            ncurve     = size(getcoef(xfdcell{ivar1}),2);
            bfdnames   = getnames(xfdcell{ivar1});
            resfdnames = bfdnames;
        else
            if any(getbasisrange(getbasis(xfdcell{ivar1})) ~= xrange)
                error('Ranges are incompatible for XFDCELL.');
            end
            if size(getcoef(xfdcell{ivar1}),2) ~= ncurve
                error('Number of curves is incompatible for XFDCELL.');
            end
        end
    end
    
    nbasmax = 0;  %  This will be the maximum number of basis functions
    
    %  check UFDCELL and extract URANGE
    %  Note:  XRANGE and URANGE are required to be identical 
    %         in this version
    
    if nforce > 0
        for ivar1=1:nvar
            urange1 = getbasisrange(getbasis(ufdcell{ivar1,1}));
            for k=1:nforce
                if ~isa_fd(ufdcell{ivar1,k})
                    error(['UFDCELL{',num2str(ivar1),',',num2str(k), ...
                            '} is not a functional data object.']);
                end  
                if ivar1 == 1 && k == 1
                    %  check that URANGE is equal to XRANGE
                    if any(urange1 ~= xrange)
                        error('XRANGE and URANGE are not identical');
                    end
                else
                    urangek = getbasisrange(getbasis(ufdcell{ivar1,k}));
                    if any(urange1 ~= urangek)
                        error('Ranges are incompatible for UFDCELL.');
                    end
                end
            end
        end
        
        %  check AWTCELL and extract the max. no. basis fns.
        
        for ivar1=1:nvar
            for k=1:nforce
                afdi = getfd(awtcell{ivar1,k});
                if ~isa_fd(afdi)
                    error('AFDI is not a functional data object.');
                end  
                basisi = getbasis(afdi);
                if any(getbasisrange(basisi) ~= urange)
                    error('Ranges are incompatible for AWTCELL.');
                end
                nbasi   = getnbasis(basisi);
                nbasmax = max([nbasmax,nbasi]);
            end
        end       
    end
    
    %  check BWTCELL
    
    if norder > 0
        for ivar1=1:nvar
            for ivar2=1:nvar
                for j=1:norder
                    if norder == 1
                        bfdPar12 = bwtcell{ivar1,ivar2};
                    else
                        bfdPar12 = bwtcell{ivar1,ivar2,j};
                    end
                    if ~isa_fdPar(bfdPar12)
                        error(['BWTCELL{',num2str(ivar1), ', ', ...
                             num2str(ivar2), ', ',              ...
                             num2str(j),                        ...
                             '} is not a functional parameter object.']);
                    end
                    bfd12  = getfd(bfdPar12);
                    basisi = getbasis(bfd12);
                    if any(getbasisrange(basisi) ~= xrange)
                        error('Ranges are incompatible in BWTCELL.');
                    end
                    nbasi   = getnbasis(basisi);
                    nbasmax = max([nbasmax,nbasi]);
                end  
            end
        end
    end
    
    %  set up sampling values to be used in numerical integration
    %    and set up matrix of basis values.  The number of sampling
    %  NFINE is here set to a usually workable value if too small.
    
    if nfine < 5*nbasmax, nfine = 5*nbasmax;  end
    
    deltax = (xrange(2)-xrange(1))/(nfine-1);
    tx     = (xrange(1):deltax:xrange(2))';
    deltau = deltax;
    tu     = tx;
    
    %  set up  YARRAY to hold values of x functions 
    %  and their derivatives
    
    yarray = zeros([nfine,ncurve,nvar,nordp1]);
    for ivar1=1:nvar
        for j=1:nordp1
            yarray(:,:,ivar1,j) = eval_fd(tx, xfdcell{ivar1}, j-1);
        end
    end
    
    %  set up  UARRAY to hold values of u functions
    
    if nforce > 0
        uarray = zeros([nfine,ncurve,nforce]);
        for k=1:nforce
            uarray(:,:,k) = squeeze(eval_fd(tu, ufdcell{ivar1,k}));
        end
    end
    
    %  set up array YPROD to hold mean of products of values in YARRAY
    
    mmat  = m2ij(nvar,nordp1);
    yprod = zeros([nfine,nvar,nordp1,nvar,nordp1]);
    for m1=1:nvar*nordp1
        i1 = mmat(m1,1);
        j1 = mmat(m1,2);
        for m2=1:m1;
            i2 = mmat(m2,1);
            j2 = mmat(m2,2);
            if ncurve == 1
                yprodval = squeeze(yarray(:,1,i1,j1)).* ...
                           squeeze(yarray(:,1,i2,j2));
            else
                yprodval = mean(squeeze(yarray(:,:,i1,j1)).* ...
                    squeeze(yarray(:,:,i2,j2)),2);
            end
            yprod(:,i1,j1,i2,j2) = yprodval;
            yprod(:,i2,j2,i1,j1) = yprodval;
        end
    end
    
    %  set up array YUPROD to hold mean of u-variables u times 
    %    x functions and their derivatives
    
    onesncurve = ones(1,ncurve);
    if nforce > 0
        yuprod = zeros(nfine, nvar, nforce, nordp1);
        for k=1:nforce
            for i1=1:nvar
                for j1=1:nordp1
                    if ncurve == 1
                        yuprodval = yarray(:,1,i1,j1).*uarray(:,1,k);
                    else
                        yuprodval = ...
                            mean(squeeze(yarray(:,:,i1,j1).* ...
                                         uarray(:,:,k)),2);
                    end
                    yuprod(:,i1,k,j1) = yuprodval;
                end
            end
        end
    end
    
    clear yarray
    
    %  set up array UPROD to hold mean of products of u-variables u 
    
    if nforce > 0
        uprod = zeros(nfine, nforce, nforce);
        for k=1:nforce
            for ju=1:k
                if ncurve == 1
                    uprodval = uarray(:,1,k).*uarray(:,1,ju);
                else
                    uprodval = mean(squeeze(uarray(:,:,k).* ...
                        uarray(:,:,ju)),2);           
                end
                uprod(:,k,ju) = uprodval;
                uprod(:,ju,k) = uprodval;
            end
            clear uarray
        end
    end
    
    %  set up an index array and some arrays of 1's
    
    mmat  = m2ij(nvar,norder);  %  order varies inside variables
    onesn = ones(nfine,1);
    
    %  set up array to hold coefficients for basis expansions
    
    if nforce > 0
        aarray = zeros(nfine,nforce);  
    else
        aarray = [];
    end
    
    if norder > 0
        barray = zeros(nfine,nvar,norder);  
    else
        barray = [];
    end
    
    
    %  -------  beginning of loop through variables  ------------
    
    for ivar1=1:nvar
        
        %  get number of coefficients to be estimated for this equation
        
        % loop through u-variables
        neqns  = 0;
        for k = 1:nforce
            afdPark = awtcell{ivar1,k};
            if getestimate(afdPark) 
                neqns = neqns + getnbasis(getbasis(getfd(afdPark)));
            end
        end
        % loop through x functions and their derivatives
        for m2=1:nvar*norder
            i2 = mmat(m2,1);
            j2 = mmat(m2,2);
            if norder == 1
                bfdParij = bwtcell{ivar1,i2};
            else
                bfdParij = bwtcell{ivar1,i2,j2};
            end
            if getestimate(bfdParij)
                neqns = neqns + getnbasis(getbasis(getfd(bfdParij)));
            end
        end
        if neqns < 1
            error('Number of equations to solve is not positive.');
        end
        
        %  set up coefficient array and right side array for 
        %     linear equation
        
        cmat   = zeros(neqns, neqns);
        dmat   = zeros(neqns, 1);
        
        %  evaluate default weight functions for this variable
        
        for k=1:nforce
            afdPark = awtcell{ivar1,k};
            aarray(:,k) = eval_fd(tu, getfd(afdPark));
        end
        
        for i=1:nvar
            for j=1:norder
                if norder == 1
                    bfdParij = bwtcell{ivar1,i};
                else
                    bfdParij = bwtcell{ivar1,i,j};
                end
                barray(:,i,j) = squeeze(eval_fd(tx, getfd(bfdParij)));
            end
        end
        
        %  loop through equations, 
        %    corresponding to rows for CMAT and DMAT
        
        %  loop through equations for u-variables
        
        mi12 = 0;
        for k1 = 1:nforce
            afdPark1  = awtcell{ivar1,k1};
            if getestimate(afdPark1)
                abasisk1    = getbasis(getfd(afdPark1));
                abasismatk1 = getbasismatrix(tu, abasisk1);
                mi11 = mi12 + 1;
                mi12 = mi12 + getnbasis(abasisk1);
                indexi1 = mi11:mi12;
                %  DMAT entry for u-variable
                weightk1 = -yuprod(:,ivar1,k1,nordp1);
                dmat(indexi1) = trapzmat(abasismatk1, onesn, deltax, ...
                                         weightk1);
                %  add terms corresponding to x-derivative weights
                %  that are not estimated
                for m=1:nvar*norder;
                    i = mmat(m,1);
                    j = mmat(m,2);
                    if norder == 1
                        bfdParij   = bwtcell{ivar1,i};
                    else
                        bfdParij   = bwtcell{ivar1,i,j};
                    end
                    if ~getestimate(bfdParij)
                        weightij = -yuprod(:,ivar1,k1,j);
                        dmat(indexi1) = dmat(indexi1) + ...
                          trapzmat(abasismatk1, barray(:,ivar1,j), deltax, ...
                                   weightij);
                    end
                end
                %  loop through weight functions to be estimated,
                %    corresponding to columns for CMAT
                %  begin with u-variable -- u-variable pairs
                mi22 = 0;
                for k2=1:nforce
                    afdPark2   = awtcell{ivar1,k2};
                    %  if estimation required, modify coefficient matrix
                    if getestimate(afdPark2)
                        abasisk2    = getbasis(getfd(afdPark2));
                        abasismatk2 = getbasismatrix(tu, abasisk2);
                        weightk2    = uprod(:,k1,k2);
                        Cprod       = ...
                            trapzmat(abasismatk1, abasismatk2, deltau, ...
                                     weightk2);
                        mi21 = mi22 + 1;
                        mi22 = mi22 + getnbasis(abasisk2);
                        indexi2 = mi21:mi22;
                        cmat(indexi1,indexi2) = Cprod;
                    end
                end
                %  remaining columns: 
                %    loop through u-variable -- x-derivative pairs
                mij22 = mi22;
                for m2=1:nvar*norder;
                    i2 = mmat(m2,1);
                    j2 = mmat(m2,2);
                    if norder == 1
                        bfdParij2 = bwtcell{ivar1,i2};
                    else
                        bfdParij2 = bwtcell{ivar1,i2,j2};
                    end
                    %  if estimation required, modify coefficient matrix
                    if getestimate(bfdParij2)
                        bbasisij2    = getbasis(getfd(bfdParij2));
                        bbasismatij2 = getbasismatrix(tx, bbasisij2);
                        weightij12   = -yuprod(:,i2,k1,j2);
                        Cprod = ...
                            trapzmat(abasismatk1, bbasismatij2, deltax, ...
                                     weightij12);
                        mij21 = mij22 + 1;
                        mij22 = mij22 + getnbasis(bbasisij2);
                        indexij2 = mij21:mij22;
                        cmat(indexi1,indexij2) = Cprod;
                    end
                end
                %  add roughness penalty matrix to diaginal entries
                lambdakl = getlambda(afdPark1);
                if lambdakl > 0.0
                    Lfdobj = getLfd(afdPark1);
                    penmat = lambdakl.*eval_penalty(abasisk1, Lfdobj);
                    cmat(indexi1,indexi1) = cmat(indexi1,indexi1) + ...
                                            penmat;
                end
            end
        end
        
        %  loop through equations for x-derivatives
        
        mij12 = mi12;
        for m1=1:nvar*norder
            i1 = mmat(m1,1);
            j1 = mmat(m1,2);
            if norder == 1
                bfdParij1 = bwtcell{ivar1,i1};
            else
                bfdParij1 = bwtcell{ivar1,i1,j1};
            end
            if getestimate(bfdParij1)
                bbasisij1    = getbasis(getfd(bfdParij1));
                bbasismatij1 = getbasismatrix(tx, bbasisij1);
                mij11 = mij12 + 1;
                mij12 = mij12 + getnbasis(bbasisij1);
                indexij1 = mij11:mij12;
                %  DMAT entry 
                weightij1 = yprod(:,i1,j1,ivar1,nordp1);
                dmat(indexij1) = ...
                    trapzmat(bbasismatij1,onesn,deltax,weightij1);
                %  add terms corresponding to forcing functions with
                %  unestimated coefficients
                for k=1:nforce
                    afdPark   = awtcell{ivar1,k};
                    if ~getestimate(afdPark)
                        weightijk = yuprod(:,ivar1,k,j1);
                        dmat(indexij1) = dmat(indexij1) + ...
                            trapzmat(bbasismatij1, aarray(:,k), deltax, ...
                                     weightijk);
                    end
                end
                %  first columns of CMAT: x-derivative -- u-variable entries
                mi22 = 0;
                for k2=1:nforce
                    afdPark2   = awtcell{ivar1,k2};
                    %  if estimation required, modify coefficient matrix
                    if getestimate(afdPark2)
                        abasisk2    = getbasis(getfd(afdPark2));
                        abasismatk2 = getbasismatrix(tx, abasisk2);
                        weightk2    = -yuprod(:,i1,k2,j1);
                        Cprod = ...
                            trapzmat(bbasismatij1, abasismatk2, deltax, ...
                                     weightk2);
                        mi21 = mi22 + 1;
                        mi22 = mi22 + getnbasis(abasisk2);
                        indexi2 = mi21:mi22;
                        cmat(indexij1,indexi2) = Cprod;
                    end
                end
                %  remaining columns: x-derivative -- x-derivative pairs
                mij22 = mi22;
                for m2=1:nvar*norder;
                    i2 = mmat(m2,1);
                    j2 = mmat(m2,2);
                    if norder == 1
                        bfdParij2 = bwtcell{ivar1,i2};
                    else
                        bfdParij2 = bwtcell{ivar1,i2,j2};
                    end
                    bbasisij2    = getbasis(getfd(bfdParij2));
                    %  if estimation required, modify coefficient matrix
                    if getestimate(bfdParij2)
                        bbasismatij2 = getbasismatrix(tx, bbasisij2);
                        weightij12   = yprod(:,i1,j1,i2,j2);
                        Cprod = ...
                            trapzmat(bbasismatij1, bbasismatij2, deltax, ...
                                     weightij12);
                        mij21 = mij22 + 1;
                        mij22 = mij22 + getnbasis(bbasisij2);
                        indexij2 = mij21:mij22;
                        cmat(indexij1,indexij2) = Cprod;
                    end
                end
                %  add roughness penalty matrix to diagonal entries
                lambdaij1 = getlambda(bfdParij1);
                if lambdaij1 > 0.0
                    Lfdobj = getLfd(bfdParij1);
                    penmat = lambdaij1.*eval_penalty(bbasisij1, Lfdobj);
                    cmat(indexij1,indexij1) = cmat(indexij1,indexij1) + ...
                                              penmat;
                end
            end
        end
        
        %  solve for coefficients of basis expansions
        
        dvec = -symsolve(cmat,dmat);
        
        %  set up u-function weight functions
        
        mi2 = 0;
        for k=1:nforce
            afdPark = awtcell{ivar1,k};
            if getestimate(afdPark)
                afdk = getfd(afdPark);
                mi1 = mi2 + 1;
                mi2 = mi2 + getnbasis(getbasis(afdk));
                indexi = mi1:mi2;
                awtcell{ivar1,k} = ...
                    putfd(afdPark, putcoef(afdk, dvec(indexi)));
            else
                awtcell{ivar1,k} = afdPark;
            end
        end
        
        %  set up X-function derivative weight functions
        
        mij2 = mi2;
        for m1=1:nvar*norder;
            i1 = mmat(m1,1);
            j1 = mmat(m1,2);
            bfdParij = bwtcell{ivar1,i1,j1};
            if getestimate(bfdParij)
                bfdij = getfd(bfdParij);
                mij1 = mij2 + 1;
                mij2 = mij2 + getnbasis(getbasis(bfdij));
                indexij = mij1:mij2;
                bwtcell{ivar1,i1,j1} = ...
                    putfd(bfdParij, putcoef(bfdij, dvec(indexij)));
            else
                bwtcell{ivar1,i1,j1} = bfdParij;
            end
        end
        
    end
    
    %  --------  end of loop through variables  -------------------
    
    %  set up residual cell RESFDCELL
    
    resfdcell     = cell(nvar);
    resfdnames{2} = 'Residual function';
    resfdnames{3} = 'Residual function value';
    
    for ivar1=1:nvar
        resbasis = getbasis(xfdcell{ivar1});
        %  initialize with highest order derivative for this variable
        resmat  = eval_fd(tx, xfdcell{ivar1}, norder);
        %  add contributions from weighted u-functions
        for k=1:nforce
            amati    = eval_fd(tu, getfd(awtcell{ivar1,k}));
            umati    = eval_fd(tu, ufdcell{ivar1,k});
            resmat   = resmat + (amati*onesncurve).*umati;
        end
        %  add contributions from weighted x-function derivatives
        for m1=1:nvar*norder;
            i1 = mmat(m1,1);
            j1 = mmat(m1,2);
            bfdij  = getfd(bwtcell{ivar1,i1,j1});
            bmatij = eval_fd(tx, bfdij)*onesncurve;
            xmatij = eval_fd(tx, xfdcell{i1}, j1-1);
            resmat = resmat + bmatij.*xmatij;
        end
        %  set up the functional data object
        resfdi = smooth_basis(tx, resmat, resbasis);
        resfdi = putnames(resfdi, resfdnames);
        resfdcell{ivar1} = resfdi;
    end
    
end

%  --------------------------------------------------------------

function mmat = m2ij(nrow,ncol)
%M2IJ sets up a NROW*NCOL by 2 matrix of row-col indices associated
%  with a number of matrix entries row-wise
%  Example:  m2ij(2,3) produces
%     1     1
%     1     2
%     1     3
%     2     1
%     2     2
%     2     3
nval = nrow*ncol;
if nval > 0
    mmat = [reshape(ones(ncol,1)*(1:nrow), nval,1), ...
            reshape((1:ncol)'*ones(1,nrow),nval,1)];
else
    mmat = [];
end

%  --------------------------------------------------------------

function XtWY = trapzmat(X,Y,delta,wt)
%TRAPZMAT integrates the products of two matrices of values
%   using the trapezoidal rule, assuming equal spacing
%  X is the first  matrix of values
%  Y is the second matrix of values
%  DELTA is the spacing between argument values (one by default)
%  WT is a vector of weights (ones by default)
%
%  XtWY is a matrix of integral estimates, number of rows equal to
%  number of col of X, number of cols equal to number of cols of Y

n = size(X,1);

if size(Y,1) ~= n
    error('X and Y do not have same number of rows.');
end

%  set default arguments

if nargin < 4, wt = ones(n,1); end
if nargin < 3, delta = 1; end

if size(wt,1) ~= n
    error('X and WT do not have same number of rows.');
end

if size(wt,2) ~= 1
    error('WT is not a column vector');
end

if delta <= 0
    error('DELTA is not a positive value.');
end

wt([1,n]) = wt([1,n])/2;
wt = wt.*delta;

X = X.*(wt*ones(1,size(X,2)));
XtWY = X'*Y;

