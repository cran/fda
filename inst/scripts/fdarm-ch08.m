%
% Ramsay, Hooker & Graves (2009)
% Functional Data Analysis with R and Matlab (Springer)
%

%  Remarks and disclaimers

%  These R commands are either those in this book, or designed to 
%  otherwise illustrate how R can be used in the analysis of functional
%  data.  
%  We do not claim to reproduce the results in the book exactly by these 
%  commands for various reasons, including:
%    -- the analyses used to produce the book may not have been
%       entirely correct, possibly due to coding and accuracy issues
%       in the functions themselves 
%    -- we may have changed our minds about how these analyses should be 
%       done since, and we want to suggest better ways
%    -- the R language changes with each release of the base system, and
%       certainly the functional data analysis functions change as well
%    -- we might choose to offer new analyses from time to time by 
%       augmenting those in the book
%    -- many illustrations in the book were produced using Matlab, which
%       inevitably can imply slightly different results and graphical
%       displays
%    -- we may have changed our minds about variable names.  For example,
%       we now prefer "yearRng" to "dayrange" for the weather data.
%    -- three of us wrote the book, and the person preparing these scripts
%       might not be the person who wrote the text
%  Moreover, we expect to augment and modify these command scripts from time
%  to time as we get new data illustrating new things, add functionality
%  to the package, or just for fun.

%
% ch. 8.  Registration: Aligning Features
%         for Samples of Curves
%

%  ----------  Registration of the Berkeley female growth data  -----------

%  Set up some strings for constructing paths to folders.
%  These strings should be modified so as to provided access
%  to the specified folders on your computer.

%  Path to the folder containing the Matlab functional data analysis
%  software

fdaMPath = 'c:/Program Files/MATLAB/R2009a/fdaM';

addpath(fdaMPath)

%  Path to the folder containing the examples

examplesPath = [fdaMPath,'/examples'];

addpath(examplesPath)

%
% Section 8.1 Amplitude and Phase Variation
%

%  Figure 8.1 in this section requires that we do
%  landmark registration first.
%  This figure is therefore plotted in Section 8.3 below.

%  Figure 8.2

%  set up a fine mesh of t-values for plotting, and define mu and sigma

tvec  = linspace(-5,5,201);
mu    = linspace(-1,1,  5);
sigma = (1:5)/3;

%  Here is function Dgauss that we use below to compute values
%  of the first derivative of a Gaussian density:


%  Data to be plotted in the left panel

DpGphase = zeros(201,5);
for i = 1:5 
    DpGphase(:,i) = DGauss(tvec+mu(i), 0, 1);
end
DpGphaseMean = mean(DpGphase,2);

%  Data to be plotted in the right panel

DpGampli = zeros(201,5);
for i = 1:5 
    DpGampli(:,i) = sigma(i)*DGauss(tvec, 0, 1);
end
DpGampliMean = mean(DpGampli,2);

%  top panel
subplot(2,1,1)
phdl=plot(tvec, DpGphase, 'b-', tvec, DpGphaseMean, 'b--', ...
          [-5,5], [0,0], 'b:');
set(phdl, 'LineWidth', 2)
axis([-5,5,-0.8,0.8])
xlabel('\fontsize{13} ') 
ylabel('\fontsize{13} ')
% bottom panel
subplot(2,1,2)
phdl=plot(tvec, DpGampli, 'b-', tvec, DpGampliMean, 'b--', ...
          [-5,5], [0,0], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} ') 
ylabel('\fontsize{13} ')
axis([-5,5,-1.2,1.2])

%  get eigenvalues for each panel

eigvecphase = svd(DpGphase).^2;
disp(eigvecphase/sum(eigvecphase));
% First three = 0.55, 0.39, and 0.05

eigvecampli = svd(DpGampli).^2;
disp(eigvecampli/sum(eigvecampli))

% First singular value = 1);
% all others are num2str-off

%
% Section 8.2 Time-Warping Functions and Registration
%

%  Figure 8.3 requires landmark registration, and is set up below

%
%  Section 8.3: Landmark registration
%

%  path to the Berkeley growth data example folder

growthPath = [examplesPath,'/growth'];

addpath(growthPath)

%  load the growth data struct object from the .mat file growth.mat

load growth

%  the dimensions of the data
ncasem = growth.ncasem;
ncasef = growth.ncasef;
nage   = growth.nage;
%  the heights of the 54 girls
heightmat = growth.hgtfmat;
%  the 31 ages of measurement
age = growth.age;

%  define the range of the ages and set up a fine mesh of ages

ageRng  = [1,18];

nfine   = 101;
agefine = linspace(ageRng(1), ageRng(2), nfine);

%  Set up functional data objects for the acceleration curves
%  and their mean.  Suffix UN means 'unregistered'.

norder = 6;
nbasis = nage + norder - 2;
growbasis = create_bspline_basis(ageRng, nbasis, norder, age);

Lfdobj    = 3;          
lambda    = 10^(-0.5);
cvecf     = zeros(nbasis, ncasef);
growfd0   = fd(cvecf, growbasis);
growfdPar = fdPar(growfd0, Lfdobj, lambda);

[Wfd, beta, hgtfhatfd] = smooth_monotone(age, heightmat, growfdPar);

accelfdUN      = deriv_fd(hgtfhatfd, 2);
accelmeanfdUN  = mean(accelfdUN);
accelmatUN     = eval_fd(agefine, accelfdUN);
accelmeanmatUN = eval_fd(agefine, accelmeanfdUN);

%  plot unregistered curves

subplot(1,1,1)
plot(accelfdUN)
xlabel('\fontsize{13} Age') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([1,18,-4,3])

%  This is a MANUAL PGS spurt identification procedure requiring
%  a mouse click at the point where the acceleration curve
%  crosses the zero axis with a negative slope during puberty.
%  Here we do this only for the first 10 children.

children = 1:10;

for icase = children
    accveci = eval_fd(agefine, accelfdUN(icase));
    plot(agefine,accveci,'b-', ageRng, [0,0], 'b--')
    xlabel('\fontsize{13} Year') 
    ylabel('\fontsize{13} Height Accel.')
    title(['Case ',num2str(icase)])
    axis([1,18,-6,4])
    pause
end

%  This is an automatic PGS spurt identification procedure.
%  A mouse click advances the plot to the next case.
%  Compute PGS mid point for landmark registration.
%  Downward crossings are computed within the limits defined
%  by INDEX.  Each of the crossings within this interval
%  are plotted.  The estimated PGS center is plotted as a vertical line.

%  The choice of range of argument values (6--18) to consider
%  for a potential mid PGS location is determined by previous
%  analyses, where they have a mean of about 12 and a s.d. of 1.

%  We compute landmarks for all 54 children

index  = 1:102;  %  wide limits
nindex = length(index);
ageval = linspace(8.5,15,nindex);
PGSctr = zeros(ncasef,1);
for icase = 1:ncasef
    accveci = eval_fd(ageval, accelfdUN(icase));
    aup     = accveci(2:nindex);
    adn     = accveci(1:(nindex-1));
    indx    = (1:102);
    indx    = indx(adn.*aup < 0 & adn > 0);
    plot(ageval(2:nindex),aup,'bo', [8,18],[0,0], 'b--')
    axis([7.9,18,-6,4])
    hold on
    for j = 1:length(indx)
        indxj = indx(j);
        aupj  = aup(indxj);
        adnj  = adn(indxj);
        agej  = ageval(indxj) + 0.1*(adnj./(adnj-aupj));
        if j == length(indx)
            PGSctr(icase) = agej;
            plot([agej,agej],[-4,4],'b-')
        else 
            plot([agej,agej],[-4,4],'b:')
        end
    end
    hold off
    title(['Case ',num2str(icase)])
    pause
end

%  We use the minimal basis function sufficient to fit 3 points
%  remember that the first coefficient is set to 0, so there
%  are three free coefficients, and the data are two boundary
%  values plus one interior knot.
%  Suffix LM means 'Landmark-registered'.

PGSctrmean = mean(PGSctr);

%  Define the basis for the function W(t).

wbasisLM = create_bspline_basis(ageRng, 4, 3, [1,PGSctrmean,18]);
WfdLM    = fd(zeros(4,1),wbasisLM);
WfdParLM = fdPar(WfdLM,1,1e-12);

%  Carry out landmark registration.

[accelfdLM, warpfdLM, wfdLM] = ...
    landmarkreg(accelfdUN, PGSctr, PGSctrmean, WfdParLM, 1);

accelmeanfdLM  = mean(accelfdLM);
accelmatLM     = eval_fd(agefine, accelfdLM);
accelmeanmatLM = eval_fd(agefine, accelmeanfdLM);

%  plot registered curves

plot(agefine, accelmatLM, 'b-', agefine, accelmeanmatLM, 'r-', ...
     [PGSctrmean,PGSctrmean], [-4,3], 'r--')
xlabel('\fontsize{13} Age') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([1,18,-4,3])

% Figure 8.1

accelmatUN10     = accelmatUN(:,children);
accelmatLM10     = accelmatLM(:,children);
accelmeanfdUN10  = mean(accelfdUN(children));
accelmeanfdLM10  = mean(accelfdLM(children));
accelmeanmatUN10 = eval_fd(agefine, accelmeanfdUN10);
accelmeanmatLM10 = eval_fd(agefine, accelmeanfdLM10);

subplot(2,1,1)
plot(agefine, accelmatUN10, 'b-', agefine, accelmeanmatUN10, 'r-', ...
     [PGSctrmean,PGSctrmean], [-4,3], 'r--')
xlabel('\fontsize{13} Age') 
ylabel('\fontsize{13} cm/yr/yr')
axis([1,18,-4,3])
subplot(2,1,2)
plot(agefine, accelmatLM10, 'b-', agefine, accelmeanmatLM10, 'r-', ...
     [PGSctrmean,PGSctrmean], [-4,3], 'r--')
xlabel('\fontsize{13} Age') 
ylabel('\fontsize{13} cm/yr/yr')
axis([1,18,-4,3])

% Figure 8.3

%  plot warping functions for cases 3 and 7

warpmatLM = eval_fd(agefine, warpfdLM);

subplot(2,2,1)
plot(agefine, accelmatUN(:,3), 'b-', ...
     [PGSctrmean,PGSctrmean], [-3,1.5], 'b--', ageRng, [0,0], 'b--')
xlabel('') 
ylabel('')
axis([1,18,-3,1.5]) 
subplot(2,2,2)
plot(agefine, warpmatLM(:,3), 'b-', ...
     [PGSctrmean,PGSctrmean], ageRng, 'b--', ...
     ageRng, ageRng, 'b--')
xlabel('') 
ylabel('')
text(PGSctrmean-0.3, warpmatLM(61,3)+0.3, 'o')
axis([1,18,1,18]) 
subplot(2,2,3)
plot(agefine, accelmatUN(:,7), 'b-', ...
     [PGSctrmean,PGSctrmean], [-3,1.5], 'b--', ageRng, [0,0], 'b--')
xlabel('') 
ylabel('')
axis([1,18,-3,1.5]) 
subplot(2,2,4)
plot(agefine, warpmatLM(:,7), 'b-', ...
     [PGSctrmean,PGSctrmean], ageRng, 'b--', ...
     ageRng, ageRng, 'b--')
xlabel('') 
ylabel('')
text(PGSctrmean-0.3, warpmatLM(61,7)+0.3, 'o')
axis([1,18,1,18]) 

%  Comparing unregistered to landmark registered curves

[MS_amp, MS_pha, RSQRLM, CLM] = ...
        AmpPhaseDecomp(accelfdUN, accelfdLM, warpfdLM, [3,18]);

disp(['Total     MS = ', num2str(MS_amp+MS_pha,2)])
disp(['Amplitude MS = ', num2str(MS_amp,2)])
disp(['Phase     MS = ', num2str(MS_pha,2)])
disp(['R-squared    = ', num2str(RSQRLM,3)])
disp(['C            = ', num2str(CLM,3)])

%
%  Section 8.4: Continuous registration
%

%  Set up a cubic spline basis for continuous registration

nwbasisCR = 15;
norderCR  =  5;
wbasisCR  = create_bspline_basis(ageRng, nwbasisCR, norderCR);
Wfd0CR    = fd(zeros(nwbasisCR,ncasef),wbasisCR);
lambdaCR  = 1;
WfdParCR  = fdPar(Wfd0CR, 1, lambdaCR);

%  carry out the registration

[accelfdCR, warpfdCR, WfdCR] = ...
              register_fd(accelmeanfdLM, accelfdLM, WfdParCR);

accelmatCR = eval_fd(agefine, accelfdCR);

subplot(1,1,1)
plot(warpfdCR)

%  plot landmark and continuously registered curves for the
%  first 10 children

accelmeanfdCR10  = mean(accelfdCR(children));
accelmeanmatCR10 = eval_fd(agefine, accelmeanfdCR10);

subplot(2,1,1)
plot(agefine, accelmatLM(:,children), 'b-', ...
     agefine, accelmeanmatLM10, 'b--', ...
     [PGSctrmean,PGSctrmean], [-3,1.5], 'b--', ageRng, [0,0], 'b--')
xlabel('\fontsize{13} Age (Years)') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([ageRng,-3,1.5])
subplot(2,1,2)
plot(agefine, accelmatCR(:,children), 'b-', ...
     agefine, accelmeanmatCR10, 'b--', ...
     [PGSctrmean,PGSctrmean], [-3,1.5], 'b--', ageRng, [0,0], 'b--')
xlabel('\fontsize{13} Age (Years)') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([ageRng,-3,1.5])

%  plot all landmark and continuously registered curves

accelmeanfdCR  = mean(accelfdCR);
accelmeanmatCR = eval_fd(agefine, accelmeanfdCR);

subplot(2,1,1)
plot(agefine, accelmatLM, 'b-', ...
     agefine, accelmeanmatLM, 'b--', ...
     [PGSctrmean,PGSctrmean], [-3,1.5], 'b--', ageRng, [0,0], 'b--')
xlabel('\fontsize{13} Age (Years)') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([ageRng,-4,3])
subplot(2,1,2)
plot(agefine, accelmatCR, 'b-', ...
     agefine, accelmeanmatCR, 'b--', ...
     [PGSctrmean,PGSctrmean], [-3,1.5], 'b--', ageRng, [0,0], 'b--')
xlabel('\fontsize{13} Age (Years)') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([ageRng,-4,3])

% Figure 8.4

subplot(1,1,1)
plot(agefine, accelmatCR(:,children), 'b-', ...
     agefine, accelmeanmatCR10, 'b--', ...
     [PGSctrmean,PGSctrmean], [-3,1.5], 'b--', ageRng, [0,0], 'b--')
xlabel('\fontsize{13} Age (Years)') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([ageRng,-3,1.5])

% Figure 8.5

phdl=plot(agefine, accelmeanmatCR, 'b-');
set(phdl, 'LineWidth', 2)
hold on
plot(agefine, accelmeanmatLM, 'b-', ...
     agefine, accelmeanmatUN, 'b--', ageRng, [0,0], 'b--')
hold off
xlabel('\fontsize{13} Age (Years)') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([ageRng,-3,1.5])

%
% Section 8.5 A Decomposition into Amplitude and Phase Sums of Squares
%

%  Comparing landmark to continuously registered curves

[MS_amp, MS_pha, RSQRCR, CCR] = ...
        AmpPhaseDecomp(accelfdLM, accelfdCR, warpfdCR, [3,18]);

disp(['Total     MS = ', num2str(MS_amp+MS_pha,2)])
disp(['Amplitude MS = ', num2str(MS_amp,2)])
disp(['Phase     MS = ', num2str(MS_pha,2)])
disp(['R-squared    = ', num2str(RSQRCR,3)])
disp(['C            = ', num2str(CCR,3)])

%
% 8.6 Registering the Chinese Handwriting Data
%

%  No code for this section

%
% 8.7 Details for Functions landmarkreg and register_fd
%

help landmarkreg
help register_fd

%
% Section 8.8 Some Things to Try
%
% (exercises for the reader)

%
% Section 8.8  More to Read
%
