function gait_pca_plot(gaitpcastr, gaitmeanfd, iharm)
gaitfinetime = linspace(0,1,101)';
harmfd  = gaitpcastr.harmfd;
meanmat = eval_fd(gaitfinetime, gaitmeanfd);
fac     = sqrt(gaitpcastr.values(iharm));
percentvar = round(100*gaitpcastr.varprop(iharm));
ivar=1;
subplot(2,1,ivar)
harmvec   = eval_fd(gaitfinetime, harmfd(iharm,ivar));
meanvec   = meanmat(:,ivar);
meanplus  = meanvec + fac.*harmvec;
meanminus = meanvec - fac.*harmvec;
phdl=plot(gaitfinetime, meanvec,   '-', [0,1], [0,0], 'r:');
set(phdl, 'LineWidth', 2)
text(gaitfinetime, meanplus,  '+')
text(gaitfinetime, meanminus, '-')
ylabel('\fontsize{12} Hip')
axis([0,1,-20,90])
title(['\fontsize{13} Harmonic ',num2str(iharm),...
    ' (Percentage of variability ',  ...
    num2str(percentvar), ')'])
ivar=2;
subplot(2,1,ivar)
harmvec   = eval_fd(gaitfinetime, harmfd(iharm,ivar));
meanvec   = meanmat(:,ivar);
meanplus  = meanvec + fac.*harmvec;
meanminus = meanvec - fac.*harmvec;
phdl=plot(gaitfinetime, meanvec,   '-', [0,1], [0,0], 'r:');
set(phdl, 'LineWidth', 2)
text(gaitfinetime, meanplus,  '+')
text(gaitfinetime, meanminus, '-')
xlabel('\fontsize{12} Normalized time')
ylabel('\fontsize{12} Knee')
axis([0,1,-20,90])
