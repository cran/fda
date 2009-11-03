% function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
%           mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
%           control = list(), intercept = 1) 
%     control = do.call('glm.control', control)
%     x       = as.matrix(x)
%     xnames  = dimnames(x)[[2L]]
%     ynames  = if is.matrix(y)) rownames(y)
%                else              names(y)
    conv    = 0;
    [nobs, nvars] = size(y);
    EMPTY = nvars == 0;
    if isempty(weights)  
        weights = ones(nobs,1);
    end
    if isempty(offset)  
        offset  = zeros(nobs,1);
    end
    variance = family.variance;
    linkinv  = family.linkinv;
%     if ~is.function(variance) || ~is.function(linkinv)) 
%         error(''family' argument seems not to be a valid family object', 
%             call. = 0)
    dev_resids  = family.dev_resids;
    aic         = family.aic;
    dmudeta     = family.dmudeta;
%     unless.null = function(x, if_null) 
%                      if isempty(x) if_null
%                      else            x
%     valideta    = unless.null(family.valideta, function(eta) 1)
%     validmu     = unless.null(family.validmu,  function(mu)  1)
%     if isempty(mustart)
%         eval(family.initialize)
%     end
%     else 
%         mukeep = mustart
%         eval(family.initialize)
%         mustart = mukeep
%     end
    if EMPTY
        %  initialization if no covariates
        eta       = zeros(nobs,1) + offset;
        if ~valideta(eta) 
            error('invalid linear predictor values in empty model');
        end
        mu        = linkinv(eta);
        if ~validmu(mu) 
            error('invalid fitted means in empty model');
        end
        dev       = sum(dev_resids(y, mu, weights));
        w         = ((weights.*dmudeta(eta).^2)./variance(mu)).^0.5;
        residuals = (y - mu)./dmudeta(eta);
        good      = ones(length(residuals),1);
        boundary  = conv == 1;
        coef      = 0;
        iter      = 0;
    else
        %  --------------------------------------------------------------
        %           normal initialization for covariates present
        %  --------------------------------------------------------------
        coefold = [];
        % initial value of link function eta
        if     ~isempty(etastart)
            eta = etastart;
        elseif ~isempty(start)
            if length(start) ~= nvars
                error('START is wrong length');
            else
                coefold = start;
                eta = offset;
                if size(x,2) == 1
                    eta = eta + x.*start;
                else
                    eta = eta + x*start;
                end
            end
        else
            eta = family.linkfun(mustart);
        end
        % initial value of expectation mu
        mu = linkinv(eta);
        if ~(validmu(mu) && valideta(eta)) 
            error('starting values not value, please specify some');
        end
        %  initial previous deviance value
        devold   = sum(dev_resids(y, mu, weights));
        boundary = conv == 0;
        %  --------------------------------------------------------------
        %                      iteration loop
        %  --------------------------------------------------------------
        for iter = 1:control.maxit
            %  indices for positive weights
            good  = weights > 0;
            %  variance values
            varmu = variance(mu);
            varmu = varmu(good);
            if any(is.na(varmu)) 
                error('NAs in V(mu)');
            end
            if any(varmu == 0)  
                error('0s in V(mu)');
            end
            %  values of derivative of mu wrt eta
            dmudeta.val = dmudeta(eta);
            if any(is_nan(dmudeta.val(good))) 
                error('NaNs in d(mu)/d(eta)');
            end
            %  indices both good and with nonzero dmudeta.val values
            good = (weights > 0) & (dmudeta.val ~= 0);
            if all(~good)
                conv = 0;
                warning(['no observations informative at iteration ', ... 
                         num2str(iter)]);
                break
            end
            %  vector of predicted values
            z = eta(good) - offset(good) + ...
                (  y(good) -     mu(good))./dmudeta.val(good);
            %  root-weights
            muvariance = variance(mu);
            w = sqrt((weights(good).*dmudeta.val(good).^2)./ ...
                      muvariance(good));
            ngoodobs = nobs - sum(~good);
            %  compute LS fit
%             fit = .Fortran('dqrls', qr = x[good, ] * w, n = ngoodobs, 
%                 p = nvars, y = w * z, ny = 1L, tol = min(1e-07, 
%                   control.epsilon/1000), coefficients = double(nvars), 
%                 residuals = double(ngoodobs), effects = double(ngoodobs), 
%                 rank = integer(1L), pivot = 1L:nvars, qraux = double(nvars), 
%                 work = double(2 * nvars), PACKAGE = 'base')
            %  check coefficients for infinite values
            if any(~is.finite(fit.coefficients))
                conv = 0;
                warning(['non-finite coefficients at iteration ', ...
                        num2str(iter)]);
                break
            end
            %  check for under-identified model
            if nobs < fit.rank 
                error(['X matrix has rank ', num2str(fit.rank), ...
                       ' but only', num2str(nobs)]);
            end
            start(fit.pivot) = fit.coefficients;
            %  update eta values
            eta = drop(x*start);
            %  update mu values
            eta = eta + offset;
            mu  = linkinv(eta);
            %  update deviance values
            dev = sum(dev_resids(y, mu, weights));
            if control.trace 
                print(['Deviance = ', num2str(dev), ...
                       'Iterations - ', num2str(iter)]);
            end
            boundary = 0;
            %  control for infinite deviance value
            if abs(dev) == inf
                if isempty(coefold) 
                  error(['no valid set of coefficients has been found:',...
                         '  please supply starting values']);
                end
                warning('step size truncated due to divergence');
                ii = 1;
                while abs(dev) == inf
                  if ii > control.maxit 
                      error('inner loop 1; cannot correct step size');
                  end
                  ii = ii + 1;
                  start = (start + coefold)./2;
                  eta = drop(x*start);
                  eta = eta + offset;
                  mu  = linkinv(eta);
                  dev = sum(dev_resids(y, mu, weights));
                end
                boundary = 1;
                if control.trace 
                    print(['Step halved: new deviance = ', num2str(dev)]);
                end
            end
            %  control for invalid values of eta and mu
            if ~(valideta(eta) && validmu(mu))
                if isempty(coefold) 
                  error(['no valid set of coefficients found: ', ...
                         'please supply starting values']);
                end
                warning('step size truncated: out of bounds');
                ii = 1;
                %  control for invalid eta or invalid mu values
                while ~(valideta(eta) && validmu(mu))
                  %  check for too many step size corrections
                  if ii > control.maxit 
                    error('inner loop 2; cannot correct step size');
                  end
                  ii    = ii + 1;
                  start = (start + coefold)./2;
                  eta   = drop(x*start);
                  eta = eta + offset;
                  mu    = linkinv(eta);
                end
                boundary = 1;
                dev      = sum(dev_resids(y, mu, weights));
                if control.trace 
                  print(['Step halved: new deviance =', num2str(dev)])
                end
            end
            %  convergence test
            if abs(dev - devold)/(0.1 + abs(dev)) < control.epsilon
                conv = 1;
                coef = start;
                break
            else
                devold  = dev;
                coefold = coef;
            end
        end
        %  --------------------------------------------------------------
        %       iteration loop completed:  wrap-up calculations
        %  --------------------------------------------------------------
        if ~conv 
            warning('glm.fit: algorithm did not converge');
        end
        if boundary 
            warning('glm.fit: algorithm stopped at boundary value');
        end
%         eps = 10 * .Machine.double.eps;
        %  warning for binomial case of 0 or 1 mu values
        if strcmp(family.family,'binomial')
            if any(mu > 1 - eps) || any(mu < eps)
                warning(['glm.fit:fittedprobabilities', ...
                         '0 or 1 occurred']);
            end
        end
        %  warning of poisson case of near zero values
        if strcmp(family.family,'poisson')
            if any(mu < eps) 
                warning('glm.fit: fitted rates numerically 0 occurred');
            end
        end
        %  control for singular design
        if fit.rank < nvars 
            error('fit.rank < nvars');
%             coef[fit.pivot][seq.int(fit.rank + 1, nvars)] = NA
        end
%         xxnames   = xnames[fit.pivot]
        residuals = (y - mu)./dmudeta(eta);
        nr        = min(sum(good), nvars);
        if nr < nvars
            Rmat = diag(nvars);
            Rmat(1:nr, 1:nvars) = fit.qr(1:nr, 1:nvars);
        else
            Rmat = fit.qr(1:nvars, 1:nvars);
        end
        Rmat(size(Rmat,1) > size(Rmat,2)) = 0;
%         names(coef)      = xnames
%         colnames(fit.qr) = xxnames
%         dimnames(Rmat)   = list(xxnames, xxnames)
    end
    %  --------------------------------------------------------------
    %               assemble list of returned values
    %  --------------------------------------------------------------
    wt = zeros(nobs,1);
    wt(good) = w.^2;
%     names(residuals) = ynames
%     names(mu)        = ynames
%     names(eta)       = ynames
%     names(wt)        = ynames
%     names(weights)   = ynames
%     names(y)         = ynames
%     if ~EMPTY 
%         names(fit.effects) = c(xxnames[seq_len(fit.rank)], rep.int('', 
%             sum(good) - fit.rank))
    if intercept 
        wtdmu = sum(weights * y)/sum(weights);
    else
        wtdmu = linkinv(offset);
    end
    nulldev   = sum(dev_resids(y, wtdmu, weights));
    n.ok      = nobs - sum(weights == 0);
    nulldf    = n.ok - as.integer(intercept);
    if EMPTY
        rank = 0;
    else
        rank = fit.rank;
    end
    resdf     = n.ok - rank;
    aic.model = aic(y, n, mu, weights, dev) + 2 * rank;
    %  set up struct object to be returned
%     list(coefficients = coef, 
%          residuals = residuals, 
%          effects   = if ~EMPTY) fit.effects, 
%          R         = if ~EMPTY) Rmat, 
%          rank      = rank, 
%          qr        = if ~EMPTY) structure(
%                          fit[c('qr', 'rank', 'qraux', 'pivot', 'tol')], 
%                          class = 'qr'), 
%          family    = family, 
%          deviance  = dev, 
%          aic       = aic.model, 
%          iter      = iter, 
%          weights   = wt, 
%          df.null   = nulldf, 
%          y         = y, 
%          converged = conv, 
%          boundary  = boundary
%          fitted.values     = mu, 
%          linear.predictors = eta, 
%          null.deviance     = nulldev, 
%          prior.weights     = weights, 
%          df.residual       = resdf, 
%         )
% end

% ___________________________________________________________________________
% 
% > glm.control
% function (epsilon = 1e-08, maxit = 25, trace = 0) 
% {
%     if ~is.numeric(epsilon) || epsilon <= 0) 
%         error('value of 'epsilon' must be > 0')
%     if ~is.numeric(maxit) || maxit <= 0) 
%         error('maximum number of iterations must be > 0')
%     list(epsilon = epsilon, maxit = maxit, trace = trace)
% end

