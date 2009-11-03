function [evalarray, basisobj] = eval_fd(evalarg, fdobj, Lfdobj)
%  EVAL_FD evaluates a functional data observation at argument
%  values EVALARG.
%
%  LFDOBJ is a functional data object defining the order m
%  HOMOGENEOUS linear differential operator of the form
%  Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) +
%          \exp[w_m(t)] D^m x(t) + ...
%
%  Arguments:
%  EVALARG ... A vector of values at which all functions are to
%              evaluated.
%  FDOBJ   ... Functional data object
%  LFDOBJ  ... A linear differential operator object
%              applied to the functions that are evaluated.
%
%  Note that the first two arguments may be interchanged.
%
%  Returns:  An array of function values corresponding
%              to the evaluation arguments in EVALARG

%  Last modified 27 September 2011

%  Check arguments

if nargin < 2
    error('Number of arguments is less than 2.');
end

%  Set default arguments

nderiv = 0;
if nargin < 3, Lfdobj = int2Lfd(0); end

%  check LFDOBJ and convert an integer to Lfd if needed.

Lfdobj = int2Lfd(Lfdobj);

%  Exchange the first two arguments if the first is an FD object
%    and the second numeric

if isnumeric(fdobj) && isa_fd(evalarg)
    temp    = fdobj;
    fdobj   = evalarg;
    evalarg = temp;
end

%  check EVALARG

if isnumeric(evalarg)
sizeevalarg = size(evalarg);
if sizeevalarg(1) > 1 && sizeevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
evalarg = evalarg(:);
elseif isstruct(evalarg)
    if ~isfield(evalarg, 'pts')
        error('Argument EVALARG does not contain a field pts.');
    end
    if ~isfield(evalarg, 'z')
        error('Argument EVALARG does not contain a field z.');
    end
else
    error('Argument evalarg is neither numeric nor a struct object.');
end

%  check FDOBJ

if ~isa_fd(fdobj)
    error('Argument FD is not a functional data object.');
end

%  Extract information about the basis

basisobj = getbasis(fdobj);
type     = getbasistype(basisobj);
rangeval = getbasisrange(basisobj);

%  check that arguments are within range

if strcmp(type, 'FEM')
    
    if nderiv ~= 0
        error('Derivative values cannot be evaluated for FEM objects.');
    end
    X = evalarg(:,1);
    Y = evalarg(:,2);
    evalarray = eval_FEM_fd(X,Y,fdobj);   
    
elseif strcmp(type, 'fdVariance')

    if nargout == 1
        evalarray = eval_fdVar_fd(evalarg, fdobj); 
    else
        [evalarray, basisobj] = eval_fdVar_fd(evalarg, fdobj);
    end

else
    
    evaldim = size(evalarg);
    temp    = reshape(evalarg,prod(evaldim),1);
    temp    = temp(~(isnan(temp)));
    EPS     = 1e-14;
    if min(temp) < rangeval(1) - EPS | max(temp) > rangeval(2) + EPS
        warning('Wid1:range', ...
            ['Values in argument EVALARG are outside of ', ...
            'permitted range, and will be ignored.']);
        disp(['Min and max args:             ', ...
            num2str([min(temp), max(temp)])])
        disp(['Min and max permitted values: ', ...
            num2str([rangeval(1), rangeval(2)])])
    end
    
    %  get maximum number of evaluation values
    
    n = evaldim(1);
    
    %  Set up coefficient array for FD
    
    coef  = getcoef(fdobj);
    coefd = size(coef);
    ndim  = length(coefd);
    if ndim <= 1
        nrep = 1;
    else
        nrep = coefd(2);
    end
    if ndim <= 2
        nvar = 1;
    else
        nvar = coefd(3);
    end
    
    %  Set up array for function values
    
    if ndim <= 2
        evalarray = zeros(n,nrep);
    else
        evalarray = zeros(n,nrep,nvar);
    end
    
    %  discard values in evalarg that are out of range
    
    if ~isempty(rangeval)
        evalarg(evalarg < rangeval(1)-1e-10) = NaN;
        evalarg(evalarg > rangeval(2)+1e-10) = NaN;
    end
    
    basismat = eval_basis(evalarg, basisobj, Lfdobj);
    
    %  evaluate the functions at arguments in EVALARG
    
    if ndim <= 2
        evalarray = full(basismat*coef);
    else
        for ivar = 1:nvar
            evalarray(:,:,ivar) = full(basismat*coef(:,:,ivar));
        end
    end
    
end

% else
% 
%     %  case of evaluation values varying from curve to curve
% 
%     for i = 1:nrep
%         evalargi = evalarg(:,i);
%         if all(isnan(evalargi))
%             error(['All values are NaN for replication ',num2str(i)]);
%         end
% 
%         index = find(~isnan(evalargi)       || ...
%                      evalargi < rangeval(1) || ...
%         evalargi > rangeval(2));
%         evalargi = evalargi(index);
%         basismat = eval_basis(evalargi, basisobj, Lfdobj);
% 
%         %  evaluate the functions at arguments in EVALARG
% 
%         if ndim == 2
%             evalarray(:,i) = NaN;
%             evalarray(index, i) = basismat*coef(:,i);
%         end
%         if ndim == 3
%             for ivar = 1:nvar
%                 evalarray(:,i,nvar) = NaN;
%                 evalarray(index,i,ivar) = basismat*coef(:,i,ivar);
%             end
%         end
%     end
% end


