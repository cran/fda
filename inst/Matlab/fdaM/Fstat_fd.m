function [F,argvals] = Fstat_fd(y,yhat,argvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fstat_fd
%
% Calculates an F-type statistic from observed and predicted responses for
% a functional linear model. This is mostly used by Fperm_fd.
%
% Arguments:
%
%   - y:  observed responses, either a vector or a functional data object. 
%
%   - yhat: predicted responses, must be of the same time as y
%
%   - argvals: time points at which to evaluate the pointwise F-statistic;
%   defaults to 101 equally spaced points. Only used if y and yhat are
%   functional data objects. 
%
% Returns:
%
%   F: the Fstatistic given as var(yhat)/mean( (y-yhat) ). If y and yhat
%   are functional data objects, this is a vector giving the F statistic at
%   each point in argvals

% Last modified May 15, 2009

    if  isnumeric(yhat), yhat = reshape(yhat,numel(yhat),1); end

    if (isnumeric(y) && ~isnumeric(yhat)) || (isa_fd(y) && ~isa_fd(yhat)) 
        error('y and yhat must both be either scalars or functional data objects.');
    end


    if isa_fd(y)
        rangeobs = getbasisrange(getbasis(y));
        rangehat = getbasisrange(getbasis(yhat));

        if any(rangeobs ~= rangehat)
            error('y and yhat do not have the same range');
        end


        if nargin < 3 || isempty(argvals)
            argvals = linspace(rangeobs(1),rangeobs(2),101);
        end

        yvec = eval_fd(argvals,y);
        yhatvec = eval_fd(argvals,yhat);

        F = var(yhatvec,[],2)./mean( (yvec-yhatvec).^2,2);
    else
        yvec = y;
        yhatvec = yhat;

        F = var(yhatvec)./mean( (yvec-yhatvec).^2 );
    end