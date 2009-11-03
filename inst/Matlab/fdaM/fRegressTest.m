%  --------------------------------------------------------------------------
%               A functional response model
%  --------------------------------------------------------------------------

N = 20;  %  No_ replications

p = 4;   %  No_ covariates (including intercept)

rangeval = [0,1];  %  argument range

%  basis for data
ndatabasis = 15;
databasis  = create_bspline_basis(rangeval, ndatabasis);

%  basis for intercept
intbasis   = create_constant_basis(rangeval);

%  define coefficients for intercept and 3 functional covariates

datasig = 1;
intcoef = randn(1,N).*datasig;
Zcoef1  = randn(ndatabasis,N).*datasig;
Zcoef2  = randn(ndatabasis,N).*datasig;

%  define corresponding functional data objects

intfd = fd(intcoef, intbasis);
Zfd1  = fd(Zcoef1, databasis);
Zfd2  = fd(Zcoef2, databasis);

%  define a scalar vector

Zvec3 = randn(N,1);

%  define basis for regression coefficient function basis

nbetabasis = 5;
betabasis = create_bspline_basis(rangeval, nbetabasis);

%  define regression coefficient function coefficients

betasig = 1;
alphacoef = randn(nbetabasis,1).*betasig;
beta1coef = randn(nbetabasis,1).*betasig;
beta2coef = randn(nbetabasis,1).*betasig;
beta3coef = randn(nbetabasis,1).*betasig;

%  define intercept and regression coefficient functions

alphafd = fd(alphacoef, betabasis);
betafd1 = fd(beta1coef, betabasis);
betafd2 = fd(beta2coef, betabasis);
betafd3 = fd(beta3coef, betabasis);

%  define a fine sequence over argument range

nfine = 501;
tfine = linspace(0,1,nfine)';

%  evaluate functions for that sequence

alphamat = repmat(eval_fd(tfine, alphafd), 1, N);
betamat1 = repmat(eval_fd(tfine, betafd1), 1, N);
betamat2 = repmat(eval_fd(tfine, betafd2), 1, N);
betavec3 = repmat(eval_fd(tfine, betafd3), 1, 1);

intmat = eval_fd(tfine,intfd);
Zmat1  = eval_fd(tfine,Zfd1);
Zmat2  = eval_fd(tfine,Zfd2);

%  define matrix of errorless response values

y0mat = intmat.*alphamat + Zmat1.*betamat1 + Zmat2.*betamat2 + ...
                          betavec3 * Zvec3';

%  define the basis, coefficients, functional data object and matrix of
%  error values

errsig  = 0.2;
errcoef = randn(ndatabasis,N).*errsig;
errfd   = fd(errcoef, databasis);
errmat  = eval_fd(tfine, errfd);

%  define the response with error values

ymat = y0mat + errmat;
yfd = smooth_basis(tfine, ymat, databasis);

%  set up the covariate and regression function Cells

ZCell    = cell(p,1);
betaCell = cell(p,1);

ZCell{1} = intfd;
ZCell{2} = Zfd1;
ZCell{3} = Zfd2;
ZCell{4} = Zvec3;

betaCell{1} = fdPar(alphafd);
betaCell{2} = fdPar(betafd1);
betaCell{3} = fdPar(betafd2);
betaCell{4} = fdPar(betafd3);

%  do the regression analysis

fRegressStr = fRegress(yfd, ZCell, betaCell);

%  plot each response curve and its fit

tplt = linspace(0,1,51)';
ymat    = eval_fd(tplt, yfd);
yhatmat = eval_fd(tplt, fRegressStr.yhat);

for i = 1:N
  phdl=plot(tplt, ymat(:,i),  'b--', tplt, yhatmat(:,i), 'b-');
  set(phdl, 'LineWidth', 2)
  title(['Curve ',num2str(i)])
  pause
end

%  plot each regression function estimate and true value

betahatCell = fRegressStr.betahat;

for j=1:4
    betahatj = eval_fd(tplt,getfd(betahatCell{j}));
    betaj    = eval_fd(tplt,getfd(betaCell{j}));
    phdl=plot(tplt, betahatj, '-', tplt, betaj, 'b--');
    set(phdl, 'LineWidth', 2)
    title(['Regression function ',num2str(j)])
    pause
end


%  ------------------------------------------------------------------------
%               A scalar response model
%  ------------------------------------------------------------------------

N = 200;  %  No_ replications

p = 4;   %  No_ covariates (including intercept)

rangeval = [0,1];  %  argument range

%  basis for data

ndatabasis = 15;
databasis  = create_bspline_basis(rangeval, ndatabasis);

%  define coefficients for intercept and 3 functional covariates

datasig = 1;
Zcoef1  = randn(ndatabasis,N).*datasig;
Zcoef2  = randn(ndatabasis,N).*datasig;

%  define corresponding functional data objects

Zfd1  = fd(Zcoef1, databasis);
Zfd2  = fd(Zcoef2, databasis);

%  define a scalar vector

Zvec3 = randn(N,1);

%  define basis for regression coefficient function basis

conbasis   = create_constant_basis(rangeval);
nbetabasis = 5;
betabasis  = create_bspline_basis(rangeval, nbetabasis);

%  define regression coefficient function coefficients

betasig   = 1;
alphacoef = randn(1)*betasig;
beta1coef = randn(nbetabasis,1).*betasig;
beta2coef = randn(nbetabasis,1).*betasig;
beta3coef = randn(1)*betasig;

%  define intercept and regression coefficient functions

alphafd = fd(alphacoef, conbasis);
betafd1 = fd(beta1coef, betabasis);
betafd2 = fd(beta2coef, betabasis);
betafd3 = fd(beta3coef, conbasis);

%  define a fine sequence over argument range

nfine = 501;
tfine = linspace(0,1,nfine)';

%  evaluate functions for that sequence

betamat1 = eval_fd(tfine, betafd1);
betamat2 = eval_fd(tfine, betafd2);

Zmat1  = eval_fd(tfine,Zfd1);
Zmat2  = eval_fd(tfine,Zfd2);

%  define matrix of errorless response values

intvec = ones(N,1);

y0vec = intvec.*alphacoef + Zmat1'*betamat1./(nfine-1) + ...
                            Zmat2'*betamat2./(nfine-1) + ...
                            Zvec3.*betacoef3;

%  define the basis, coefficients, functional data object and matrix of
%  error values

errsig  = 0.1;
errvec  = randn(N,1).*errsig;

%  define the response with error values

yvec = y0vec + errvec;

%  set up the covariate and regression function Cells

ZCell    = cell(4,1);
betaCell = cell(4,1);

ZCell{1} = intvec;
ZCell{2} = Zfd1;
ZCell{3} = Zfd2;
ZCell{4} = Zvec3;

betaCell{1} = fdPar(alphafd);
betaCell{2} = fdPar(betafd1);
betaCell{3} = fdPar(betafd2);
betaCell{4} = fdPar(betafd3);

%  do the regression analysis

fRegressStr = fRegress(yvec, ZCell, betaCell);

yhat = fRegressStr.yhat;

plot(yvec, yhat, 'o', yhat, yhat, '--');

betahatCell = fRegressStr.betahat;

tplt = linspace(0,1,51);
for j=1:4
    betahatj = eval_fd(tplt,getfd(betahatCell{j}));
    betaj    = eval_fd(tplt,getfd(betaCell{j}));
    phdl=plot(tplt, betahatj, '-', tplt, betaj, 'b--');
    set(phdl, 'LineWidth', 2)
    title(['Regression function ',num2str(j)])
    pause
end
