function SSE_CV = fRegress_CV(yvec, xfdcell, betacell)
% FREGRESS_CV computes cross-validated error sum of squares

%  Last modified 20 July 2006

N = length(yvec);
p = length(xfdcell);

SSE_CV  = 0;
for i=1:N
    clear xfdcelli;
    indexi = find((1:N) ~= i);
    xfdcelli = cell(p,1);
    for j=1:p
        xfd = xfdcell{j};
        xfdcelli{j} = xfd(indexi);
    end
    yveci = yvec(indexi);
    fRegressCelli = fRegress(yveci, xfdcelli, betacell);
    betaestcelli  = fRegressCelli{4};
    yhati = 0;
    for j=1:p
        betafdj = getfd(betaestcelli{j});
        xfdj    = xfdcell{j};
        bbasisj = getbasis(betafdj);
        rangej  = getbasisrange(bbasisj);
        nfine   = max(101, getnbasis(bbasisj)*10+1);
        tfine   = linspace(rangej(1), rangej(2), nfine)';
        delta   = tfine(2)-tfine(1);
        betavec = eval_fd(tfine, betafdj);
        xveci   = eval_fd(tfine, xfdj(i));
        yhati   = yhati + delta.*(sum(xveci.*betavec) - ...
            0.5.*( xveci(1)    *betavec(1) + ...
                   xveci(nfine)*betavec(nfine) ));
    end
    SSE_CV = SSE_CV + (yvec(i) - yhati).^2;
end



