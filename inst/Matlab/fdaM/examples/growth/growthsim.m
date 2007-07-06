
addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\growmatlab')
addpath ('c:\matlab\fdaM\examples\growth')

%  Last modified 9 January 2004

%  -----------------------------------------------------------------------
%      Simulate Berkeley Growth Data using the Jolicoeur Model
%                  as the true curve
%  -----------------------------------------------------------------------

%   set up the ages of observation for Berkeley data 

age   = [ 1:0.25:2, 3:8, 8.5:0.5:18 ]';
n     = length(age);
rng   = [1,18];
nones = ones(1,n);

%  set up the functional data object for the standard error of 
%  measurement.  The coefficients were computed empirically from 
%  analyses of the the Fels data.

stderrrng   = [0,22];
stderrbasis = create_bspline_basis(stderrrng, 7, 4);
stderrcoef  = [ ...
   0.7702    0.6961    0.4606    0.5066    0.5161    0.4700    0.5168]';
stderrfd = fd(stderrcoef, stderrbasis);

%  plot the standard error of measurement

plot(stderrfd)

agefine = linspace(1,18,101)';
stderrfine = eval_fd(agefine, stderrfd);

subplot(1,1,1)
plot(agefine, stderrfine, 'k-')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Standard error of measurement (cm)')

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/FelsStdErr.ps'

%  evaluate this standard error of measurement at the ages 
%  of observation

stderrvec = eval_fd(age, stderrfd);

%  The Jolicoeur model is:
%  h(t) = a*[1 - 1/F(t)] where
%    F(t) = [b_1*(t + e)]^c_1 + [b_2*(t + e)]^c_2 + [b_3*(t + e)]^c_3

%  Set up the mean vector and variance matrix for the coefficients
%  These were obtained empirically from analyses of the Fels data

%  the coefficients are in the order:

%          a, b_1, b_2, b_3, c_1, c_2, c_3, e

%  set up the mean vector for the coefficients
Coefmean = [ ...
 164.7077, 0.3071, 0.1106, 0.0816, 0.7285, 3.6833, 16.6654, 1.4744]';

%  set up the standard deviations for the coefficients
Coefstddev  = [ ...
    5.8892, 0.0425, 0.0078, 0.0058, 0.0591, 0.2232, 0.7424, 0.3183]';

%  set up the correlations among the coefficients
Coefcor  = [ ...
    1.0000   -0.1437   -0.0536   -0.1659    0.0008    0.1033   -0.0607   -0.1173;
   -0.1437    1.0000    0.3522    0.4710   -0.1490    0.5203   -0.2652   -0.8659;
   -0.0536    0.3522    1.0000    0.5832   -0.0948   -0.0868   -0.1022   -0.3000;
   -0.1659    0.4710    0.5832    1.0000    0.3923    0.3278   -0.7472   -0.2920;
    0.0008   -0.1490   -0.0948    0.3923    1.0000    0.1039   -0.4341    0.1949;
    0.1033    0.5203   -0.0868    0.3278    0.1039    1.0000   -0.4481   -0.4837;
   -0.0607   -0.2652   -0.1022   -0.7472   -0.4341   -0.4481    1.0000    0.1523;
   -0.1173   -0.8659   -0.3000   -0.2920    0.1949   -0.4837    0.1523    1.0000];

%  compute the variance-covariance matrix of the coefficients

Coefvar = diag(Coefstddev) * Coefcor * diag(Coefstddev);

%  compute the Choleski factor for the coefficients
%  An N by 8 coefficient matrix can be simulated by computing
%  Z*Coeffac, where Z is an N by 8 matrix of independent standard
%  normal coefficients

Coeffac =  chol(Coefvar);

%  simulate a sample of N Jolicoeur model curves

N = 10;  
Z       = randn(N,8);
Coefmat = Z*Coeffac + ones(N,1)*Coefmean';

%  evaluate the heights, velocities and accelerations

hgttrue = zeros(N,n);
veltrue = zeros(N,n);
acctrue = zeros(N,n);
for i=1:N
    coefi = Coefmat(i,:)';
    [hgti, veli, acci] = jolifn(age, coefi);
    hgttrue(i,:) = hgti';
    veltrue(i,:) = veli';
    acctrue(i,:) = acci';
end

%  add some error with the appropriate standard deviation

hgtobs = hgttrue;
for i=1:N
    errorvec = randn(1,n).*stderrvec';
    hgtobs(i,:) = hgtobs(i,:) + errorvec;
end

%  plot observed and true curves

subplot(1,1,1)
for i=1:N
    plot(age, hgtobs(i,:), 'o', age, hgttrue(i,:), 'b-')
    xlabel('\fontsize{16} Age')
    ylabel('\fontsize{16} Height (cm)')
    title(['Case ',num2str(i)])
    pause
end

%  ---------------------------------------------------------------
%         Smoothing the data using spline smoothing
%  --------------------------------------------------------------- 

%  set up the basis for smoothing the data

norder = 6;
nbasis = norder + n - 2;
basis  = create_bspline_basis(rng, nbasis, norder, age);

%  smooth the data

lambda   = 1e-1;
Lfd      = 4;
hgtfdPar = fdPar(basis, Lfd, lambda);

wtvec  = 1.0./stderrvec.^2;
wtvec  = wtvec./mean(wtvec);

[hgtfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_basis(hgtobs', age, hgtfdPar, wtvec);

%  compute mean squared errors

accest = eval_fd(age, hgtfd, 2);
MSE    = zeros(n,1);
for i=1:N
    [hgttruei, veltruei, acctruei] = jolifn(age, Coefmat(i,:)');
    MSE = MSE + (accest(:,i) - acctruei).^2./N;
end

[df, gcv, MSE([11,19])']

%  plot the estimated and the true accelerations

subplot(1,1,1)
agefine = linspace(1,18,101)';
accest  = eval_fd(agefine, hgtfd, 2);
MSE = zeros(101,1);
for i=1:N
    coefi = Coefmat(i,:)';
    [hgttruei, veltruei, acctruei] = jolifn(agefine, coefi);
    accesti = accest(:,i);
    MSE = MSE + (accesti - acctruei).^2./N;
    plot(agefine, accesti, '-', agefine, acctruei, '--')
    xlabel('\fontsize{16} Age')
    ylabel('\fontsize{16} Acceleration (cm/sec^2)')
    title(['Case ',num2str(i)])
    pause
end

plot(agefine, MSE)

%   Smooth using a range of smoothing parameter values

%  set up coefficient matrix for model

N = 1000;
Z = randn(N,8);
Coefmat = Z*Coeffac + ones(N,1)*Coefmean';

%  generate the true curves and random data

hgttrue = zeros(N,n);
veltrue = zeros(N,n);
acctrue = zeros(N,n);
hgtobs  = zeros(N,n);
for i=1:N
    coefi = Coefmat(i,:)';
    [hgti, veli, acci] = jolifn(age, coefi);
    hgttrue(i,:) = hgti';
    veltrue(i,:) = veli';
    acctrue(i,:) = acci';
    errorvec     = randn(1,n).*stderrvec';
    hgtobs(i,:)  = hgttrue(i,:) + errorvec;
end

%  plot a sample of true accelerations

accfine = zeros(20,101);
agefine = linspace(1,18,101)';
for i=1:20
    coefi = Coefmat(i,:)';
    [hgti, veli, acci] = jolifn(agefine, coefi);
    accfine(i,:) = acci';
end

subplot(1,1,1)
plot(agefine, accfine, '-', [1,18], [0,0], 'r:')
xlabel('\fontsize{19} Age')
ylabel('\fontsize{19} Acceleration (cm/sec^2)')
axis([1,18,-5,2.2])

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/JoliAccel.ps'

%  set up the range of log lambda values

loglam  = (-2.5:0.25:0.5)';
nlam    = length(loglam);
dfsave  = zeros(nlam,1);
gcvsave = zeros(nlam,1);
MSEsave = zeros(nlam,n);
Varsave = zeros(nlam,n);
Biasave = zeros(nlam,n);

%  loop through the log lambda values

for ilam=1:nlam
    lambda = 10^loglam(ilam);
    hgtfdPar = fdPar(basis, Lfd, lambda);
    [hgtfd, df, gcv] = smooth_basis(age, hgtobs', hgtfdPar, wtvec);
    accest = eval_fd(age, hgtfd, 2);
    dfsave(ilam)    = df;
    gcvsave(ilam)   = gcv;
    MSEsave(ilam,:) = mean((accest'-acctrue).^2);
    Varsave(ilam,:) = var(accest'-acctrue);
    Biasave(ilam,:) = mean(accest'-acctrue);
end

%  display degrees of freedom values

[loglam,dfsave,gcvsave]

%  plot GCV and root mean squared errors for ages 8, 12, and 16

index = [11,19,27];  % indices of ages 8, 12, and 16
subplot(2,1,1)
plot(loglam, gcvsave, '-', [-1,-1], [400,600], 'r:')
ylabel('\fontsize{16} GCV')
subplot(2,1,2)
plot(loglam, sqrt(MSEsave(:,11)), 'g-', ...
     loglam, sqrt(MSEsave(:,19)), 'b-', ...
     loglam, sqrt(MSEsave(:,27)), 'r-', ...
     [-1,-1], [0,1], 'r:')
xlabel('\fontsize{16} log_{10} \lambda')
ylabel('\fontsize{16} Accel. RMSE')
legend('\fontsize{16} 8','12','16')

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/GCVMSEplot.ps'

%  plot squared bias, variance, and MSE for three ages

subplot(1,1,1)
for i=index
    plot(loglam, MSEsave(:,i),    '-', ...
         loglam, Biasave(:,i).^2, '-', ...
         loglam, Varsave(:,i),    '-')
    xlabel('\fontsize{16} log_{10} \lambda')
    title(['\fontsize{16} Results for age ',num2str(age(i))])
    legend('\fontsize{12} Mean squared error', ...
           'Squared Bias', 'Sampling Variance')
    pause
end

%  plot root mean squared error over all ages for best lambda

ilam = 7

subplot(3,1,1)
plot(age, sqrt(MSEsave(ilam,:)), '-')
ylabel('\fontsize{16} RMSE')
axis([1,18,0,2])

%  plot bias

subplot(3,1,2)
plot(age, Biasave(ilam,:), '-', [1,18], [0,0], 'r:')
ylabel('\fontsize{16} Bias')

%  plot standard error of estimate

subplot(3,1,3)
plot(age, sqrt(Varsave(ilam,:)), '-')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Std. Err.')
axis([1,18,0,1.1])

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/MSEBiaVar.ps'

%  ---------------------------------------------------------------
%         Display pointwise confidence limits for acceleration
%                      using spline smoothing
%  --------------------------------------------------------------- 

%  set up the basis for smoothing the data

norder = 6;
nbasis = norder + n - 2;
basis  = create_bspline_basis(rng, nbasis, norder, age);

%  smooth the data

lambda   = 1e-1;
Lfd      = 4;
hgtfdPar = fdPar(basis, Lfd, lambda);

wtvec  = 1.0./stderrvec.^2;
wtvec  = wtvec./mean(wtvec);

[hgttrue, veltrue, acctrue] = jolifn(age, Coefmean);

[hgtfd, df, gcv, coef, SSE, penmat, y2cMap] = ...
    smooth_basis(age, hgttrue, h, hgtfdPar, wtvec);

%  compute mean squared errors

accbasis  = eval_basis(agefine, basis, 2);
Lmat      = accbasis*y2cMap;
Varacc    = diag(Lmat * diag(stderrvec) * Lmat');
stderracc = sqrt(Varacc);

[hgtfine, velfine, accfine] = jolifn(agefine, Coefmean);

subplot(1,1,1)
plot(agefine, accfine, 'k-', ...
     agefine, accfine + 2.*stderracc, 'k--', ...
     agefine, accfine - 2.*stderracc, 'k--', ...
     [1,18], [0,0], 'k:')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Height acceleration')

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/AccStdErr.ps'

subplot(1,1,1)
plot(age, Lmat( 54,:), 'k-', ...
     age, Lmat(  1,:), 'k--', ...
     age, Lmat(101,:), 'k-.', ...
     [1,18], [0,0], 'k:')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Acceleration weight')
axis([1,18,-1.3,1.3])
legend('\fontsize{12} Age 10', 'Age 1', 'Age 18')

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/AccWeight.ps'

%  ---------------------------------------------------------------
%      Display pointwise confidence limits for acceleration
%                 using monotone smoothing
%  --------------------------------------------------------------- 

%  set up the basis for smoothing the data

norder   = 6;
nbasis   = norder + n - 2;
basis    = create_bspline_basis(rng, nbasis, norder, age);

Lfd      = int2Lfd(2);
lambda   = 1e-10;
hgtfdPar = fdPar(basis, Lfd, lambda);

Wfd0   = smooth_basis(age, acctrue./veltrue, hgtfdPar);


%  set up monotone smooth

conv    = 1e-3;            %  convergence criterion
iterlim = 20;              %  max. no. iterations
active  = 1:nbasis;        %  indices of coefficients to be optimized
dbglev  = 0;               %  level of output per iteration (full here)
wtvec   = 1./stderrvec.^2; %  weight vector for smoothing
wtvec   = wtvec./mean(wtvec);
zmat    = ones(n,1);

%  set hgtfdPar

Lfd      = 3;               %  penalize curvature of acceleration
lambda   = 10^(-1);         %  smoothing parameter
hgtfdPar = fdPar(Wfd0, Lfd, lambda);

%  smooth the data

N = 100;
hgtest  = zeros(N,n);
velest  = zeros(N,n);
accest  = zeros(N,n);
coefest = zeros(N,nbasis);
betaest = zeros(N,2);
hwait  = waitbar(0,'Please wait...');
for i=1:N
    errorvec = randn(n,1).*stderrvec;
    hgtobs   = hgttrue + errorvec;
    [Wfd, beta] = ...
        smooth_monotone(age, hgtobs, hgtfdPar, zmat, wtvec, ...
                        conv, iterlim, active, dbglev);
    hgtest(i,:) = beta(1) + beta(2).*eval_mon(age, Wfd, 0)';
    velest(i,:) =           beta(2).*eval_mon(age, Wfd, 1)';
    accest(i,:) =           beta(2).*eval_mon(age, Wfd, 2)';
    coefest(i,:) = getcoef(Wfd)';
    betaest(i,:) = beta';
    waitbar(i/N,hwait)
end

save growthsim

%  compute mean squared errors

accres = accest - ones(N,1)*acctrue';
accmse = mean(accres.^2)';
stderrmonacc = sqrt(accmse);

accbasis  = eval_basis(age, basis, 2);
Lmat      = accbasis*y2cMap;
Varacc    = diag(Lmat * diag(stderrvec) * Lmat');
stderracc = sqrt(Varacc);

plot(age, stderrmonacc, '-', age, stderracc, '--')

subplot(1,1,1)
plot(agefine, accfine, 'k-', ...
     age,     acctrue + 2.*stderracc,    'k--', ...
     age,     acctrue - 2.*stderracc,    'k--', ...
     age,     acctrue + 2.*stderrmonacc, 'k--', ...
     age,     acctrue - 2.*stderrmonacc, 'k--', ...
     [1,18], [0,0], 'k:')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Height acceleration')
axis([1,18,-6,3])
for i=2:n-1
    line([age(i),age(i)], ...
         [acctrue(i) - 2.*stderrmonacc(i), 
          acctrue(i) + 2.*stderrmonacc(i)])
end

print -dps2 'c:/MyFiles/fdabook1/revision/figs.dir/AccStdMonErr.ps'

%  set up confidence intervals

%  smooth the errorless curve
hgtfine = jolifn(agefine, Coefmean');

lambda = 1e-1;
[Wfd, beta, Fstr, iternum, iterhist, y2cMap] = ...
    smooth_monotone(age, hgttrue, hgtfdPar, zmat, wtvec, ...
                    conv, iterlim, active, dbglev);
Wfd0 = Wfd;
CVar   = (y2cMap * diag(stderrvec.^2) * y2cMap');
coef   = getcoef(Wfd);
coefse = sqrt(diag(CVar));

% [coef, coefse]

basismat = eval_basis(age, getbasis(Wfd));
WVar     = basismat*CVar*basismat';

Wstderr = sqrt(diag(WVar));

Wvec    = eval_fd(age, Wfd);
plot(age, Wvec, '-', ...
     age, Wvec - 2.*Wstderr, 'g--', ...
     age, Wvec + 2.*Wstderr, 'g--', ...
     [1,18],[0,0],'r:')
axis([1,18,-3,2])
     
Wfdest = fd(coefest', basis);
Wfdmat = eval_fd(age, Wfdest);
Wfdlo = zeros(n, 1);
Wfdhi = zeros(n, 1);
for i=1:n
    Wfdveci = sort(Wfdmat(i,:))';
    Wfdlo(i) = Wfdveci(floor(N*0.05)+1);
    Wfdhi(i) = Wfdveci(floor(N*0.95)+1);
end

plot(age, Wvec, '-', ...
     age, Wvec - 2.*Wstderr, 'g--', ...
     age, Wvec + 2.*Wstderr, 'g--', ...
     age, Wfdlo,             'm--', ...
     age, Wfdhi,             'm--', ...
     [1,18], [0,0], 'r:')
axis([1,18,-3,2])

%  ---------------------------------------------------------------
%         Simulate a sample of N Jolicoeur model curves,
%      smoothing the data using monotone spline smoothing
%  --------------------------------------------------------------- 

%  simulate a sample of N Jolicoeur model curves

N       = 10;
Z       = randn(N,8);
Coefmat = Z*Coeffac + ones(N,1)*Coefmean';

%  evaluate the heights, velocities and accelerations

hgttrue = zeros(N,n);
veltrue = zeros(N,n);
acctrue = zeros(N,n);
for i=1:N
    coefi = Coefmat(i,:)';
    [hgti, veli, acci] = jolifn(age, coefi);
    hgttrue(i,:) = hgti';
    veltrue(i,:) = veli';
    acctrue(i,:) = acci';
end

%  add some error with the appropriate standard deviation

hgtobs = hgttrue;
for i=1:N
    errorvec = randn(1,n).*stderrvec';
    hgtobs(i,:) = hgtobs(i,:) + errorvec;
end

%  set up the basis for smoothing the data

norder = 5;
nbasis = norder + n - 2;
basis  = create_bspline_basis(rng, nbasis, norder, age);
Wfd0   = fd(zeros(nbasis,1),basis);
zmat   = ones(n,1);

%  smooth the data

Lfd     = 3;               %  penalize curvature of acceleration
conv    = 1e-3;            %  convergence criterion
iterlim = 10;              %  max. no. iterations
active  = 1:nbasis;        %  indices of coefficients to be optimized
dbglev  = 0;               %  level of output per iteration (full here)
lambda  = 10^(-1);         %  smoothing parameter

% [hgtfd, df, gcv] = smooth_basis(hgtobs', age, basis, ...
%                                 stderrvec.^2, Lfd, lambda);

accest = zeros(N,n);
for i=1:N
    disp(i)
    [Wfd, beta] = ...
        monotone(age, hgtobs(i,:)', stderrvec.^2, ...
        Wfd0, zmat, Lfd, lambda, ...
        conv, iterlim, active, dbglev);
    accest(i,:) = beta(2).*eval_mon(age, Wfd, 2)';
end

%  compute mean squared errors
MSE    = zeros(n,1);
for i=1:N
    [hgttruei, veltruei, acctruei] = jolifn(age, Coefmat(i,:)');
    MSE = MSE + (accest(i,:)' - acctruei).^2./N;
end
MSE([11,19])'


%  plot the estimated and the true accelerations

subplot(1,1,1)
agefine = linspace(1,18,101)';
for i=1:N
    coefi = Coefmat(i,:)';
    [hgttruei, veltruei, acctruei] = jolifn(age, coefi);
    accesti = accest(i,:)';
    plot(age, accesti, '-', age, acctruei, '--')
    xlabel('\fontsize{16} Age')
    ylabel('\fontsize{16} Acceleration (cm/sec^2)')
    title(['Case ',num2str(i)])
    pause
end

%   Smooth using a range of smoothing parameter values

%  set up coefficient matrix for model

N = 100;
Z = randn(N,8);
Coefmat = Z*Coeffac + ones(N,1)*Coefmean';

%  generate the true curves and random data

hgttrue = zeros(N,n);
veltrue = zeros(N,n);
acctrue = zeros(N,n);
hgtobs  = zeros(N,n);
for i=1:N
    coefi = Coefmat(i,:)';
    [hgti, veli, acci] = jolifn(age, coefi);
    hgttrue(i,:) = hgti';
    veltrue(i,:) = veli';
    acctrue(i,:) = acci';
    errorvec     = randn(1,n).*stderrvec';
    hgtobs(i,:)  = hgttrue(i,:) + errorvec;
end

%  plot true accelerations

plot(age, acctrue)
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Acceleration (cm/sec^2)')

%  set up the range of log lambda values

loglam  = (-3:0.25:0)';
nlam    = length(loglam);
MSEmonsave = zeros(nlam,n);
Varmonsave = zeros(nlam,n);
Biamonsave = zeros(nlam,n);

%  loop through the log lambda values

for ilam=1:nlam
    lambda = 10^loglam(ilam);
    %     [hgtfd, df, gcv] = smooth_basis(hgtobs', age, basis, ...
    %                                     stderrvec.^2, Lfd, lambda);
    accest = zeros(N,n);
    for i=1:N
        fprintf('\n%3.f %5.f', [ilam,i])
        [Wfd, beta] = ...
            monotone(age, hgtobs(i,:)', stderrvec.^2, ...
            Wfd0, zmat, Lfd, lambda, ...
            conv, iterlim, active, dbglev);
        accest(i,:) = beta(2).*eval_mon(age, Wfd, 2)';
    end
    MSEmonsave(ilam,:) = mean((accest-acctrue).^2);
    Varmonsave(ilam,:) = var(accest-acctrue);
    Biamonsave(ilam,:) = mean(accest-acctrue);
end

%  plot GCV and root mean squared errors for ages 8, 12, and 16

index = [11,19,27];  % indices of ages 8, 12, and 16
subplot(1,1,1)
plot(loglam, sqrt(MSEmonsave(:,index)), '-', ...
     [-1.5,-1.5], [0,1], 'b:')
xlabel('\fontsize{12} log_{10} \lambda')
ylabel('\fontsize{12} Accel. RMSE')
legend('8','12','16')

%  plot squared bias, variance, and MSE for three ages

subplot(1,1,1)
for i=index
    plot(loglam, MSEmonsave(:,i), '-', ...
         loglam, Biamonsave(:,i).^2, '-', ...
         loglam, Varmonsave(:,i), '-')
    xlabel('\fontsize{16} log_{10} \lambda')
    title(['\fontsize{16} Results for age ',num2str(age(i))])
    legend('\fontsize{12} Mean squared error', ...
           'Squared Bias', 'Sampling Variance')
    pause
end

%  plot root mean squared error over all ages for best lambda

ilam = 7;
ilammon = 4;

plot(age, sqrt(MSEmonsave(ilammon,:)), '-', ...
     age, sqrt(MSEsave(ilam,:)), '--')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Root mean squared error')

%  plot bias

plot(age, Biamonsave(ilammon,:), '-', ...
     age, Biasave(ilam,:), '-')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Bias')

%  plot standard error of estimate

plot(age, sqrt(Varmonsave(ilammon,:)), '-', ...
     age, sqrt(Varsave(ilam,:)), '-')
xlabel('\fontsize{16} Age')
ylabel('\fontsize{16} Std. Error of Estimate')

