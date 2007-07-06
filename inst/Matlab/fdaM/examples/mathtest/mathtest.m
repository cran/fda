addpath('c:\matlab\fdaM')  %  add path to FDA functions 

%  ----------------------------------------------------------------
%  ----------------------  enter ACT male data  -------------------
%  ----------------------------------------------------------------

nit =   60;
nex = 2115;

%  Read in the data.  The first record contains the key,
%    and subsequent records the responses for each examinee

fid = fopen('actm.txt','rt');
ACTmtest = reshape(fscanf(fid, '%s'), [nit,nex])';

%  set up limits on quadrature points

thetamax = floor(20*(1 + log10(nex))/3)/10;
thetamin = -thetamax;
thetarng = [thetamin, thetamax];

%  Smooth the data using the TestGraf algorithm that
%  estimates item response functions by kernel smoothing
%  The objective is to get an initial estimate of the
%  item response functions

nbin = 51;  %  number of bin boundaries
thetavec = linspace(thetamin, thetamax, nbin)';  %  bin boundaries

%  bin the data

Pbin0  = FirstStep(ACTmtest, thetavec);

%  set up smoothing matrix used in TestGraf

Smat = smoothmat(thetavec, nex);

%  smooth the binned data

Pbin = Smat*Pbin0;

%  compute log odds of smoothed data

Wbin  = log(Pbin./(1-Pbin)); 

%   set up B-spline basis  

nbasis = 11;
norder = 6;
wbasis = create_bspline_basis(thetarng, nbasis, norder);

%  set up a functional data object for initial W

Wfd0  = data2fd(Wbin, thetavec, wbasis); 

%  compute asymptotic standard error

Wvec  = eval_fd(Wfd0,0);
DWvec = eval_fd(Wfd0,1);
Pvec  = exp(Wvec)./(1+exp(Wvec));
serr  = 1/sqrt(sum(DWvec.^2.*Pvec.*(1-Pvec)));

%  compute number of quadrature points

nq = floor((thetamax - thetamin)/serr);  % nq = 26

%  set up equally-spaced trait values theta_q and thetafine

thetaq    = linspace(thetamin, thetamax, nq)';  %  quadrature points
thetafine = linspace(thetamin, thetamax, 101)'; %  plotting   points

%  least squares weights for B-spline test functions

norder = 4;
nbasis = nq + norder - 2;
wgtq   = gausswgt(thetaq, nbasis, norder);

% Matrix phimat has basis function values for each thetaq value

phimat = getbasismatrix(thetaq, wbasis);

%  set up linear differential operator

Ldata = zeros(101,3);
theta = linspace(thetamin,thetamax,101)';
Ldata(:,3) = (exp(theta)-1)./(exp(theta)+1);
wfd = data2fd(Ldata,theta,wbasis);
Lfdobj = Lfd(3,fd2cell(wfd));

% Matrix K is the penalty matrix for penalizing roughness;

Kmat = eval_penalty(wbasis, Lfdobj);  %Penalize size of 2nd derivative

%  ------------------------------------------------------------------
%                         The EM Algorithm
%  ------------------------------------------------------------------

%  Function FirstStep computes the proportions of examinees
%    associated with each value of thetaq that pass items

Pbin  = FirstStep(ACTmtest, thetaq);
Wbin  = log(Pbin./(1-Pbin)); % matrix of logit transforms of probabilities
Wfd0  = data2fd(Wbin, thetaq, wbasis); % functions corresponding to these
coef0 = getcoef(Wfd0);         %  coefficients for the basis expansion
W0    = eval_fd(Wfd0, thetaq); % evaluate functions at thetaq values
P0    = 1./(1+exp(-W0));    %  compute corresponding probabilities
P     = P0; Q = 1-P;        % starting values for probabilities
coef  = coef0;              % starting values for coefficients
Wfd   = Wfd0;               % starting functions w

%  ------------  initialize the EM algorithm  --------

lambda   = 1e0;              % smoothing parameter
penmat   = lambda.*Kmat;      % penalty matrix times smoothing parameter
iter     = 0;                 % initialize iteration number
convtest = 1e-3;              % convergence criterion
itermax  = 100;                % maximum number of iterations
F        = 1e10;              % initialize new and old function values
Fold     = F + 2*convtest;

%  ---- iterate through EM algorithm  calling functions Estep and Mstep --------

while Fold - F > convtest & iter < itermax
   if iter == 0, disp('It.    -log L   penalty    F     Fold - F'); end
   iter = iter + 1;
   Fold = F;

  %  --------  E step  ------------

  [N, CN, CP, L, CL] = Estep(ACTmtest, P, wgtq);
   
  %  compute penalized negative sum of marginal likelihoods
  logL = sum(log(L));
  pen  = sum(diag(coef' * penmat * coef));
  F = -logL + pen;
  %  print out F, which is being minimized
  fprintf('%g  ', [iter, -logL, pen, F, Fold - F]);  fprintf('\n');
   
  %   -------  Mstep  -------------
  
  coef = Mstep(CN, N, P, coef, phimat, penmat);
  % update Q by n matrix of probabilities 
  P = 1./(1+exp(-phimat * coef));  
  Q = 1 - P;
  
end

Wfd = putcoef(Wfd, coef);

%  --------------------  end of EM algorithm ---------------------

%  set results for fine mesh of theta values

Wfdfine = eval_fd(thetafine, Wfd);
Pfine   = exp(Wfdfine)./(1 + exp(Wfdfine));
Qfine   = 1 - Pfine;

%  compute arc length along manifold 

DWfdmat = eval_fd(thetaq, Wfd, 1);
DPmat   = P.*Q.*DWfdmat;
DPsqr   = sum((DPmat').^2);
DPnorm  = sqrt(DPsqr)';

DWfdfine   = eval_fd(thetafine, Wfd, 1);
DPfine     = Pfine.*Qfine.*DWfdfine;
DPnormfine = sqrt(sum((DPfine').^2))';

Smax  = (thetaq(2) - thetaq(1)).*(sum(DPnorm) - ...
         0.5.*(DPnorm(1) + DPnorm(nq)));
Msfd  = data2fd(log(DPnorm), thetaq, wbasis);
Svec  = monfn(thetaq,Msfd);
Svec  = Smax.*Svec./max(Svec);
Sfine = Smax.*monfn(thetafine,Msfd)./monfn(thetamax,Msfd);

DsPmat  = DPmat./(DPnorm*ones(1,nit));
DsPfine = DPfine./(DPnormfine*ones(1,nit));

%  estimate density function for trait scores

sumCP = sum(CP);
pdf   = wgtq'.*sumCP;
pdf   = pdf./sum(pdf');

%  express W as a function of arc length

Wmat    = eval_fd(thetaq, Wfd);
Wsbasis = create_bspline_basis([0,max(Svec)], 25);
Wsfd    = data2fd(Wmat, Svec, Wsbasis);

%  --------------------------------------------------------
%                     Display results
%  --------------------------------------------------------

%  plot "data" and fit for all items

subplot(1,1,1)
itemindex = 1:nit;
for j = itemindex
  plot(thetaq, (CN(j,:)./N)', '.', ...
       thetaq, P(:,j), 'b-', thetaq, Pbin(:,j), 'go')
  axis([thetamin, thetamax, 0, 1])
  xlabel('\fontsize{16} \theta')
  ylabel('\fontsize{16} P(\theta)')
  title(['\fontsize{16} Item ',num2str(j)])
  pause
end

itemindex = [1, 9, 59];
m = 0;
for j = itemindex
    m = m + 1;
    subplot(1,3,m)
    hdl = plot(thetaq, P(:,j), 'k-');
    set(hdl, 'LineWidth', 2)
    axis([thetamin, thetamax, 0, 1])
    axis('square')
    xlabel('\fontsize{16} \theta')
    title(['\fontsize{16} Item ',num2str(j)])  
    if j == 1
        ylabel('\fontsize{16} Prob. correct')
    end
end

print -dps2 'c:\MyFiles\fdabook1\revision\figs.dir\threeirfs.ps'

% plot log odds functions

plot(Wfd)
xlabel('\fontsize{16} \theta')
ylabel('\fontsize{16} W(\theta)')

%  plot arc length

subplot(1,1,1)
plot(thetafine, Sfine, 'k-')
axis([-3,3,0,9])
xlabel('\fontsize{16} \theta')
ylabel('\fontsize{16} Arc Length s')

print -dps2 'figs.dir\actarclength.ps'

%  plot squared slopes against arc length one by one

itemindex = 1:nit;
for j = itemindex
  plot(Svec, DsPmat(:,j).^2)
  axis([0, 10, 0, .25])
  xlabel('\fontsize{16} Arc length s')
  ylabel('\fontsize{16} Squared Slope')
  title(['\fontsize{16} Item ',num2str(j)])
  pause
end

%  plot all squared slopes against arc length

itemindex=1:nit;
plot(Sfine, DsPfine(:,itemindex).^2, 'k-')
axis([0, 8, 0, .2])
xlabel('\fontsize{16} Arc length s')
ylabel('\fontsize{16} Squared Slope')

print -dps2 'figs.dir\actdiscrim.ps'

%  plot estimated information functions for items 1, 9 and 59

itemindex=[1,9,59];
info = DsPmat(:,itemindex).^2./(P(:,itemindex).*(1-P(:,itemindex)));
plot(Svec, info)
axis([0, 9.5, 0, 1])
xlabel('\fontsize{16} Arc length s')
ylabel('\fontsize{16} Squared Slope')

%  plot results for item 56

subplot(1,2,1)
plot(Sfine, Pfine(:,56))
axis([0,Smax,0,1])
axis('square')
xlabel('\fontsize{12} Arc length s')
ylabel('\fontsize{12} Probability of Success')

subplot(1,2,2)
plot(Sfine, DsPfine(:,56).^2)
axis([0,Smax,0,0.2])
axis('square')
xlabel('\fontsize{12} Arc length s')
ylabel('\fontsize{12} Squared Slope')

%  plot conditional likelihoods in CL for selected examinees

subplot(1,1,1)
index = 1:10;
for i = index 
  CLmax = max(CL(i,:));
  plot(thetaq,CL(i,:)./CLmax, '-')
  xlabel('\fontsize{16} \theta')
  ylabel('\fontsize{16} Conditional Likelihood')
  title(['\fontsize{16} Examinee ',num2str(i)])
  pause
end

%  plot density function for ability

subplot(1,1,1)
plot(thetaq, pdf)
xlabel('\fontsize{16} \theta')
ylabel('\fontsize{16} Density')
title(['\fontsize{16} Probability Density Function for ', ...
       'Trait Score Values'])

%  plot arc length against theta

subplot(1,1,1)
plot(thetaq, Svec)
xlabel('\fontsize{16} \theta')
ylabel('\fontsize{16} Arc Length')

print -dps2 'arclength.ps'

% plot density function for arclength values

plot(Svec, pdf'./DPnorm)
xlabel('\fontsize{16} Arc Length')
ylabel('\fontsize{16} Density')
title(['\fontsize{16} Probability Density Function for ', ...
       'Arc Length Values'])

%  plot probability manifold for three items

subplot(1,1,1) 
index3 = [1,9,59]; % indices of three items to be plotted
%  plot curve
P1 = P(:,index3(1));
P2 = P(:,index3(2));
P3 = P(:,index3(3));
plot3(P1, P2, P3, 'o-')
axis([0,1,0,1,0,1])
xlabel(['\fontsize{16} Item ',num2str(index3(1))])
ylabel(['\fontsize{16} Item ',num2str(index3(2))])
zlabel(['\fontsize{16} Item ',num2str(index3(3))])
%  enable rotation by mouse dragging
grid on
rotate3d

print -dps2 'figs.dir\actmanifold.ps'

%  Usual plot of same IRF's both as P and as W functions

W1 = log(P1./(1-P1));
W2 = log(P2./(1-P2));
W3 = log(P3./(1-P3));

subplot(1,2,1)
plot(thetaq, P1, '-', thetaq, P2, '.-', thetaq, P3, 'o-')
axis([thetamin, thetamax, 0, 1])
axis('square')
text(0, P1((nq)/2)+.05, ' 1')
text(0, P2((nq)/2)+.1, ' 9')
text(0, P3((nq)/2)+.1, '59')
xlabel('\fontsize{16} Ability \theta')
title('\fontsize{16} Probability of Success P(\theta)')

subplot(1,2,2)
plot(thetaq, W1, '-', thetaq, W2, '.-', thetaq, W3, 'o-')
axis('square')
xlabel('\fontsize{16} Ability \theta')
title('\fontsize{16} Log Odds-Ratio W(\theta)')

%  principal components analysis of W functions

nharm   = 4;
Wpcastr = pca(Wfd, nharm);
Wpcastr = varmx_pca(Wpcastr);

Wharmmat = eval_fd(thetaq, Wpcastr.harmfd);

Wfdmean  = mean(Wfd);
Wmean    = eval_fd(thetaq, Wfdmean);
Wmeanmat = Wmean*ones(1,4);

Pmean   = exp(Wmean)./(1 + exp(Wmean));

Wconst = ones(nq,1)*[2, 1, .5, .5];
Wmatp = Wmeanmat +  Wconst.*Wharmmat;
Wmatm = Wmeanmat -  Wconst.*Wharmmat;
Pharp = exp(Wmatp)./(1+exp(Wmatp));
Pharm = exp(Wmatm)./(1+exp(Wmatm));

titlestr = ['  I: 25%'; ' II: 15%'; 'III: 24%'; ' IV: 35%'];
for j=1:nharm
  subplot(2,2,j)
  plot(thetaq, Pmean, '--')
  axis([thetamin, thetamax, 0, 1])
  text(thetaq-.1, Pharp(:,j), '+')
  text(thetaq-.1, Pharm(:,j), '-')
  title(['\fontsize{12} Harmonic ',titlestr(j,:)])
end

print -dps2 'figs.dir\actPCA'

%  plot components one by one

subplot(1,1,1)
Wplot_pca(Wpcastr);






