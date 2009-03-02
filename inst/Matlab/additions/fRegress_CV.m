function [SSE_CV,errfd] = fRegress_CV(yfdPar, xfdcell, betacell,cvobs)
% FREGRESS_CV computes cross-validated error sum of squares

%  Last modified 20 July 2006

if strcmp(class(yfdPar),'double')

    yvec = yfdPar;
    N = length(yvec);
    p = length(xfdcell);

    if nargin<4,
        cvobs = 1:N;
        N2 = N;
    else
        N = length(cvobs);
    end

    SSE_CV  = 0;
    for m=1:N2
        i = cvobs(m);
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
elseif isa_fdPar(yfdPar) || isa_fd(yfdPar)
    if isa_fdPar(yfdPar), yfd = getfd(yfdPar); else yfd = yfdPar; end

    ycoef = getcoef(yfd);
    N = size(ycoef,2);
    if nargin < 4,
        cvobs = 1:N;
        N2= N;
    else
        N2 = length(cvobs);
    end;

    p = length(xfdcell);

    SSE_CV = 0;
    errcoefs = [];

    for m = 1:N2
        i =  cvobs(m);
        if m == 2
            print(strcat('Estimated Computing time =  ',num2str(round(N*elapsed.time)),'seconds.'))
        end

        tic
        clear xfdcelli;
        indexi = find((1:N) ~= i);
        xfdcelli = cell(p,1);
        for j=1:p
            xfd = xfdcell{j};
            xfdcelli{j} = xfd(indexi);
        end

        yfdi = yfd(indexi);
        fRegressCelli = fRegress(yfdi, xfdcelli, betacell);

        yhat = 0;                        % Now we predict
        for k = 1:p
            xfdk = xfdcell{k};
            yhat = yhat + xfdk(i).*getfd(fRegressCelli{4}{k});
        end
        err = yfd(i) - yhat;

        errcoefs = [errcoefs,getcoef(err)];

        SSE_CV = SSE_CV + inprod(err,err);
        elapsed.time = toc;
    end
    errfd = fd(errcoefs,getbasis(err));
    %    names(errfd$fdnames)[[3]] = "Xval Errors"    %% NEED TO CHECK THIS


else
    %  YFDOBJ is neither functional nor multivariate
    error('YFDOBJ is neither functional nor multivariate.');
end


