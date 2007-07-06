function Wplot_pca(pcastr, nx, pointplot, harm, expand, cycle)
%  WPLOT_PCA  Plots the harmonics for a principal components analysis
%    of W functions for test.
%  Arguments:
%  PCASTR    ... Struct object returned by PCA_FD.
%  NX        ... Number of argument values for plotting. Default = 128.
%  POINTPLOT ... If pointplot=1, then the harmonics are plotted with
%                +'s and -'s on either side of the mean function.
%                Otherwise lines are used.  Default is 1.
%  HARM      ... If harm = 0 (the default) then all the computed harmonics
%                are plotted.   Otherwise those in HARM are plotted.
%  EXPAND    ... If expand =0 then effect of +/- 2 standard deviations of
%                each harmonic are given.
%                Otherwise the factor expand is used.
%  CYCLE     ... If cycle=T and there are 2 variables then a cycle plot
%                will be drawn.  If the number of variables is anything else,
%                CYCLE will be ignored.

%  last modified 14 Sept 2000

  %  set up default argument values

  if nargin < 6
    cycle = 0;
  end
  if nargin < 5
    expand = 0;
  end
  if nargin < 4
    harm = 0;
  end
  if nargin < 3
    pointplot = 0;
  end
  if nargin < 2
    nx = 101;
  end

  harmfd    = pcastr.harmfd;
  basisobj  = getbasis(harmfd);
  rangex    = getbasisrange(basisobj);
  fdnames   = getnames(harmfd);
  x         = linspace(rangex(1), rangex(2), nx);
  fdmat     = eval_fd(x, harmfd);
  meanmat   = squeeze(eval_fd(x, pcastr.meanfd));
  Wmeanmat  = exp(meanmat);
  Wmeanmat  = Wmeanmat./(1+Wmeanmat);
  dimfd     = size(fdmat);
  nharm     = dimfd(2);
  if harm == 0
    harm = (1:nharm);
  end
  if length(dimfd) == 2
    for iharm = harm
      if expand == 0
        fac = sqrt(pcastr.eigvals(iharm));
      else
        fac = expand;
      end
      vecharm = fdmat(:,iharm);
      percentvar = round(100 * pcastr.varprop(iharm));
      meanplus   = meanmat+fac.*vecharm;
      meanminus  = meanmat-fac.*vecharm;
      Wmeanplus  = exp(meanplus);
      Wmeanplus  = Wmeanplus./(1+Wmeanplus);
      Wmeanminus = exp(meanminus);
      Wmeanminus = Wmeanminus./(1+Wmeanminus);
      plottop = 1;
      plotbot = 0;
      plot(x, Wmeanmat,   '-')
      text(x, Wmeanplus,  '+')
      text(x, Wmeanminus, '-')
      xlabel(fdnames{1});
      ylabel(fdnames{3});
      axis([rangex(1), rangex(2), plotbot, plottop])
      title(['PCA function ', num2str(iharm), ...
             ' (Percentage of variability ', num2str(percentvar), ')'])
      disp('Press any key to continue')
      pause;
    end
  end

