function plot_pca(pcastr, nx, harm, expand, cycle)
%  PLOT_PCA  Plots the harmonics for a principal components analysis.
%  Arguments:
%  PCASTR    ... Struct object returned by PCA_FD.
%  NX        ... Number of argument values for plotting. Default = 101.
%  HARM      ... If harm = 0 (the default) then all the computed harmonics
%                are plotted.   Otherwise those in HARM are plotted.
%  EXPAND    ... If expand =0 then effect of +/- 2 standard deviations of
%                each harmonic are given.
%                Otherwise the factor expand is used.
%  CYCLE     ... If cycle=T and there are 2 variables then a cycle plot
%                will be drawn.  If the number of variables is anything else,
%                CYCLE will be ignored.

%  last modified 25 July 2006

  %  set up default argument values

  if nargin < 5
    cycle = 0;
  end
  if nargin < 4
    expand = 0;
  end
  if nargin < 3
    harm = 0;
  end
  if nargin < 2
    nx = 101;
  end

  harmfd  = pcastr.harmfd;
  basis   = getbasis(harmfd);
  rangex  = getbasisrange(basis);
  fdnames = getnames(harmfd);
  x       = linspace(rangex(1), rangex(2), nx);
  fdmat   = eval_fd(harmfd, x);
  meanmat = squeeze(eval_fd(pcastr.meanfd, x));
  dimfd   = size(fdmat);
  nharm   = dimfd(2);
  if harm == 0
    harm = (1:nharm);
  end
  if length(dimfd) == 2
    %  plotting for univariate functions
    for iharm = harm
      if expand == 0
        fac = sqrt(pcastr.values(iharm));
      else
        fac = expand;
      end
      vecharm = fdmat(:,iharm);
      percentvar = round(100 * pcastr.varprop(iharm));
      meanplus  = meanmat+fac.*vecharm;
      meanminus = meanmat-fac.*vecharm;
      plottop = max([meanplus;meanminus]);
      plotbot = min([meanplus;meanminus]);
      plot(x, meanmat,   '-')
      text(x, meanplus,  '+')
      text(x, meanminus, '-')
      xlabel(fdnames{1});
      ylabel(fdnames{3});
      axis([rangex(1), rangex(2), plotbot, plottop])
      title(['PCA function ', num2str(iharm), ...
             ' (Percentage of variability ', num2str(percentvar), ')'])
      disp('Press any key to continue')
      pause;
    end
  else
    %  plotting for multivariate functions
    if cycle && dimfd(3) == 2
      %  cycle plotting for  bivariate functions
      for iharm = harm
        if expand == 0
          fac = 2 * sqrt(pcastr.values(iharm));
        else
          fac = expand;
        end
        vecharm = fdmat(:,iharm);
        percentvar = round(100 * pcastr.varprop(iharm));
        meanplus  = meanmat+(fac.*vecharm) * ones(1,2);
        meanminus = meanmat-(fac.*vecharm) * ones(1,2);
        plottop = max(max([meanplus, meanminus, meanmat]));
        plotbot = min(min([meanplus, meanminus, meanmat]));
        plot(meanmat(:,1),   meanmat(:,2),   '-')
        text(meanplus(:,1),  meanplus(:,2),  '+')
        text(meanminus(:,1), meanminus(:,2), '-')
        xlabel(fdnames{1});
        ylabel(fdnames{3});
        axis([plotbot, plottop, plotbot, plottop])
        title(['PCA function ', num2str(iharm), ...
               ' (Percentage of variability ', num2str(percentvar), ')'])
        disp('Press any key to continue')
        pause;
      end
    else
      for iharm = harm
        if expand == 0
          fac = sqrt(pcastr.values(iharm));
        else
          fac = expand;
        end
        nvar = dimfd(3);
        subplot(nvar,1,1)
        for jvar = 1:nvar
          percentvar = round(100 * pcastr.varprop(iharm));
          vecharm = fdmat(:,iharm,jvar);
          meanplus  = meanmat(:,jvar)+fac.*vecharm;
          meanminus = meanmat(:,jvar)-fac.*vecharm;
          plottop = max([meanplus;meanminus]);
          plotbot = min([meanplus;meanminus]);
          subplot(nvar,1,jvar)
          plot(x, meanmat(:,jvar), '-')
          text(x, meanplus,  '+')
          text(x, meanminus, '-')
          xlabel(fdnames{1});
          ylabel(fdnames{3});
          axis([rangex(1), rangex(2), plotbot, plottop])
          title(['Harmonic ', num2str(iharm), ...
             ' Variable ',    num2str(jvar),  ...
             ' (Percentage of variability ', num2str(percentvar), ')'])
        end
        disp('Press any key to continue')
        pause;
      end
    end
  end

