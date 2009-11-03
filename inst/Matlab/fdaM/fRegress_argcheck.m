function [yfdPar, xfdcell, betacell, wt, rangeval] = ...
               fRegress_argcheck(yfdPar, xfdcell, betacell, wt)
%  FREGRESS_ARGCHECK checks the first four arguments for the functions
%  for function regression, including FREGRESS.

%  --------------------  Check classes of arguments  --------------------

%  check YFDPAR and compute sample size N

if isa_fd(yfdPar)
    yfdPar = fdPar(yfdPar);
end

if ~(isa_fdPar(yfdPar) || strcmp(class(yfdPar), 'double'))
    error('Argument YFDPAR is not of class fd, fdPar or double.');
end

if isa_fdPar(yfdPar)
    yfd   = getfd(yfdPar);
    ycoef = getcoef(yfd);
    N     = size(ycoef,2);
end

if isnumeric(yfdPar)
    N = length(yfdPar);
end

%  default weights

if nargin < 4 || isempty(wt),  wt = ones(N,1);  end

%  Check that xfdcell is a cell object

if isa_fd(xfdcell) || isnumeric(xfdcell)
    xfdcell = {xfdcell};
end

if ~iscell(xfdcell)
    error('Argument XFDCELL is not a cell object.');
end

%  get number of independent variables 
    
p = length(xfdcell);

%  Check BETACELL

if isa_fd(betacell)
    betacell = {betacell};
end
  
if ~iscell(betacell)
    error('Argument BETACELL is not a cell object.');
end

if length(betacell) ~= p
    error(['Number of regression coefficients does not match', ...
           ' number of independent variables.']);
end

%  check that the regression is functional, and extract the range

if isa_fdPar(yfdPar)
    rangeval = getbasisrange(getbasis(getfd(yfdPar)));
else
    allscalar = 1;
    for j=1:p
        if isa_fd(xfdcell{j})
            rangeval = getbasisrange(getbasis(xfdcell{j}));
            allscalar = 0;
            break
        end
    end
    if allscalar
        error(['The dependent variable and all the independent ', ...
               'variables are scalar.']);
    end
end

%  --------------------  check contents of arguments  -------------------

%  If the object is a vector of length N,
%  it is converted to a functional data object with a
%  constant basis

onebasis = create_constant_basis(rangeval);
onesfd   = fd(1,onebasis);

berror = 0;
xerror = 0;
for j=1:p
    
    %  XFDCELL:
    
    xfdj = xfdcell{j};
    if isa_fd(xfdj)
        xcoef = getcoef(xfdj);
        if length(size(xcoef)) > 2
            error(['Covariate ',j,' is not univariate.']);
        end
        %  check size of coefficient array
        Nj = size(xcoef,2);
        if Nj ~= N
            disp(['Incorrect number of replications in XFDLIST ', ...
                'for covariate ', num2str(j)])
            xerror = 1;
        end
    end
    if isnumeric(xfdj)
        xvecj = xfdj;
        Zdimj = size(xvecj);
        if Zdimj(1) ~= N && Zdimj(1) ~= 1
            disp(['Vector in XFDLIST{', num2str(j), ...
                '} has wrong length.'])
            xerror = 1;
        end
        if Zdimj(2) ~= 1
            disp(['Matrix in XFDLIST{', num2str(j), ...
                '} has more than one column.'])
            xerror = 1;
        end
        if size(xvecj,1) == 1
            xvecj = xvecj*ones(N,1);
        end
        xfdj = fd(xvecj', onebasis);
        xfdcell{j} = xfdj;
    end
    if ~(isa_fd(xfdj) || isnumeric(xfdj))
        disp(['XFDCELL{', num2str(j), ...
            '} is neither an FD object nor numeric.'])
        xerror = 1;
    else
        %  check range
        xrngj = getbasisrange(getbasis(xfdj));
        if xrngj(1) ~= rangeval(1) || xrngj(2) ~= rangeval(2)
            disp(['Range incompatibility in XFDCELL{', num2str(j), '}.'])
            xerror = 1;
        end
    end
    
    %  BETACELL:
    
    betafdParj = betacell{j};
    if isa_fd(betafdParj) || isa_basis(betafdParj)
        betafdParj  = fdPar(betafdParj);
        betacell{j} = betafdParj;
    end
    if ~isa_fdPar(betafdParj)
        disp(['BETACELL{',num2str(j),'} is not a FDPAR object.']);
        berror = 1;
    else
        %  check range
        brngj = getbasisrange(getbasis(getfd(betafdParj)));
        if brngj(1) ~= rangeval(1) || brngj(2) ~= rangeval(2)
            disp(['Range incompatibility in BETACELL{', num2str(j), '}.'])
            berror = 1;
        end
    end
end

if xerror || berror
    error('An error has been found in either XFDCELL or BETACELL.'); 
end

%  WT:

if isempty(wt)
    wt = ones(N,1);
end
   
if length(wt) ~= N 
    error('Number of weights not equal to N.');
end

if any(wt < 0)     
    error('Negative weights found.');
end
