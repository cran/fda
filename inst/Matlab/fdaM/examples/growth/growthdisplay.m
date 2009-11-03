function growthdisplay(age, hgtmat, hgtfd, stderr, y2cMap, gender)
%GROWTHDISPLAY displays:
%         height curves with the raw data
%         residuals from height curves
%         velocity     curves with 95% confidence bands
%         acceleration curves with 95% confidence bands

%  Arguments are:
%  AGE     ...  The ages at which measurements are taken
%  HGTMAT  ...  A matrix of raw data for height
%  HGTFD   ...  A functional data object for height
%  STDERR  ...  A vector of standard error measurement values for
%               each age
%  Y2CMAP  ...  A matrix produced by SMOOTH_BASIS that maps
%               a data vector into the corresponding coefficient
%               vector
%  GENDER  ...  A string, either 'male' or 'female'

%  Last modified 31 October 2003

%  get number of ages and number of cases
[nage,ncase] = size(hgtmat);
%  get range of ages and set up a fine set of values for plotting
agelo = min(age);
agehi = max(age);
agefn = linspace(agelo, agehi, 101)';  %  ages for plotting
rng   = [agelo, agehi];

%  evaluate heights for each age
hgtfit = eval_fd(age,     hgtfd);
%  evaluate height, velocity, and acceleration for each 
%    plotting value of age
hgthat = eval_fd(agefn, hgtfd);
velhat = eval_fd(agefn, hgtfd, 1);
acchat = eval_fd(agefn, hgtfd, 2);

%  set up the standard errors of velocity and acceleration

hgtbasis  = getbasis(hgtfd);

velbasis  = eval_basis(agefn, hgtbasis, 2);
velLmat   = velbasis*y2cMap;
stderrvel = sqrt(diag(velLmat * diag(stderr) * velLmat'));

accbasis  = eval_basis(agefn, hgtbasis, 2);
accLmat   = accbasis*y2cMap;
stderracc = sqrt(diag(accLmat * diag(stderr) * accLmat'));

%  do the plotting for each case, pausing after each plot

for i = 1:ncase
    %  plot height and raw data
    subplot(2,2,1)
    plot(age, hgtmat(:,i), 'o', ...
         agefn,  hgthat(:,i), '-')
    axis([agelo, agehi,60,200]);
    xlabel('Years');  
    title(['Height for ',gender,' ',num2str(i)])
    %  plot residuals
    subplot(2,2,2)
    resi = hgtmat(:,i) - hgtfit(:,i);
    ind  = find(resi >= -.7 & resi <= .7);
    plot(age(ind), resi(ind), 'o-', rng, [0,0], '--')
    axis([agelo,agehi,-.7,.7])
    xlabel('Years');  title('Residuals')
    %  plot velocity curve with confidence limits
    subplot(2,2,3)
    ind = find(velhat(:,i) >= 0 & velhat(:,i) <= 20);
    plot(agefn(ind), velhat(ind,i), '-', ...
        agefn(ind), velhat(ind,i) + 2.*stderrvel(ind), 'b:', ...
        agefn(ind), velhat(ind,i) - 2.*stderrvel(ind), 'b:')
    axis([agelo,agehi,0,20])
    xlabel('Years');  title('Velocity')
    %  plot acceleration curve with confidence limits
    subplot(2,2,4)
    ind = find(acchat(:,i) >= -6 & acchat(:,i) <= 6);
    plot(agefn(ind), acchat(ind,i), '-', ...
        agefn(ind), acchat(ind,i) + 2.*stderracc(ind), 'b:', ...
        agefn(ind), acchat(ind,i) - 2.*stderracc(ind), 'b:', ...
        rng, [0,0], 'r--')
    axis([agelo,agehi,-6,6]),
    xlabel('Years');  title('Acceleration')
    pause
end
