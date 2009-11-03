function [SSE_CV, errfd] = ...
              fRegress_CV(yfdPar, xfdcell, betacell, wt, CVobs)
% FREGRESS_CV computes cross-validated error sum of squares
% for scalar or functional responses. NOTE: ordinary and
% generalized cross validation scores are now returned by fRegress
% when scalar responses are used.

%  Last modified 30 September 2009 by Jim Ramsay

if nargin < 4, wt = [];  end

[yfdPar, xfdcell, betacell, wt] = ...
               fRegress_argcheck(yfdPar, xfdcell, betacell, wt);
          
p = length(xfdcell);
N = size(getcoef(xfdcell{1}),2);

if nargin < 5, CVobs = 1:N;  end
M = length(CVobs);

if isnumeric(yfdPar)  
    %  Dependent variable is scalar
    yvec   = yfdPar;
    SSE_CV = 0;
    errvec = zeros(M,1);
    for m=1:M
        i        = CVobs(m);
        indexi   = find((1:N) ~= i);
        if ~isempty(wt)
            wti  = wt(indexi);
        else
            wti  = [];
        end
        xfdcelli = cell(p,1);
        for j=1:p
            xfd = xfdcell{j};
            xfdcelli{j} = xfd(indexi);
        end
        yveci         = yvec(indexi);
        fRegressStri = fRegress(yveci, xfdcelli, betacell, wti);
        betaestcelli  = fRegressStri.betahat;
        yhati = 0;
        for j=1:p
            betafdj = getfd(betaestcelli{j});
            xfdj    = xfdcell{j};
            bbasisj = getbasis(betafdj);
            rangej  = getbasisrange(bbasisj);
            nfine   = max(501, getnbasis(bbasisj)*10+1);
            tfine   = linspace(rangej(1), rangej(2), nfine)';
            delta   = tfine(2)-tfine(1);
            betavec = eval_fd(tfine, betafdj);
            xveci   = eval_fd(tfine, xfdj(i));
            yhati   = yhati + delta.*(sum(xveci.*betavec) - ...
                0.5.*( xveci(1)    *betavec(1) + ...
                       xveci(nfine)*betavec(nfine) ));
        end
        errvec(i) = yvec(i) - yhati;
        SSE_CV    = SSE_CV + errvec(i).^2;
    end
    errfd = errvec;
else
    %  Dependent variable is functional
    yfd      = getfd(yfdPar);
    SSE_CV   = 0;
    errcoefs = [];
    for m = 1:M
        i = CVobs(m);
        indexi = find((1:N) ~= i);
        if ~isempty(wt)
            wti  = wt(indexi);
        else
            wti  = [];
        end
        xfdcelli = cell(p,1);
        for j=1:p
            xfd = xfdcell{j};
            xfdcelli{j} = xfd(indexi);
        end
        yfdi = yfd(indexi);        
        fRegresscelli = fRegress(yfdi,xfdcelli,betacell,wti);
        betaestcelli = fRegresscelli.betahat;
        yhatfdi = 0;
        for j=1:p
            betafdParj = betaestcelli{j};
            betafdj    = getfd(betafdParj);
            xfdj       = xfdcell{j};
            xfdij      = xfdj(i);
            yhatfdi    = yhatfdi + xfdij.*betafdj;
        end
        errfdi   = yfd(i) - yhatfdi;
        SSE_CV   = SSE_CV + inprod(errfdi,errfdi);
        errcoefs = [errcoefs, getcoef(errfdi)];
    end
    if nargout > 1
        errfd = fd(errcoefs,getbasis(errfdi));
    end
end



