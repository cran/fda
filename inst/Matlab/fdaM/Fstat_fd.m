function [F, argvals] = Fstat_fd(y,yhat,argvals)
%  FSTAT_FD computes an F-ratio scalar or functional data object

if nargin < 3, argvals = [];  end

% check arguments

if class(y) ~= class(yhat)
    error('Y and YHAT are not of the same class.');
end

if ~(isa_fd(y) || isa_double(y))
    error('Y and YHAT are neither numeric nor FD objects.');
end

%  branch depending on the class of Y

if isa_fd(y)
    rangeobs = getbasisrange(getbasis(y));
    rangehat = getbasisrange(getbasis(yhat));
    if ~all(rangeobs == rangehat)
        error('Y and YHAT do not have the same range');
    end
    if isempty(argvals)
        argvals = linspace(rangeobs(1),rangeobs(2),101);
    end
    yvec    = eval_fd(argvals,y);
    yhatvec = eval_fd(argvals,yhat);
    F = var(yhatvec,0,2)./mean((yvec-yhatvec).^2,2);
else
    yvec    = y;
    yhatvec = yhat;
    F       = var(yhatvec)/mean((yvec-yhatvec).^2);
end


