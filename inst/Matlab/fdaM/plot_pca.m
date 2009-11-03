function plot_pca_fd(pcastr, matplt, harm, expand, cycle, colorwrd)
%  PLOT_PCA  Plots the harmonics for a principal components analysis.
%  Arguments:
%  PCASTR    ... Struct object returned by PCA_FD.
%  MATPLT    ... If MATPLT = 0 (the default), then the harmonics
%                are displayed one by one with a keystroke to page
%                between them.  Otherwise they are displayed together.
%  HARM      ... If HARM = 0 (the default) then all the computed harmonics
%                are plotted.   Otherwise those in HARM are plotted.
%  EXPAND    ... If expand =0 then effect of +/- 2 standard deviations of
%                each harmonic are given.
%                Otherwise the factor expand is used.
%  CYCLE     ... If cycle=T and there are 2 variables then a cycle plot
%                will be drawn.  If number of variables is anything else,
%                CYCLE will be ignored.
%   COLORWRD ... If 0, plus and minus harmonic values are plotted as
%                symbols + and - resp., but otherwise, as colors
%                green and red.

%  This is identical to plot_pca_fd and is included for 
%  compatibility with the corresponding R function.

%  last modified 16 September 2009

%  set up default argument values

if nargin < 6
    colorwrd = 1;
end
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
    matplt = 0;
end

harmfd    = pcastr.harmfd;
basisobj  = getbasis(harmfd);
nbasis    = getnbasis(basisobj);
nx        = max([101, 10*nbasis+1]);
rangex    = getbasisrange(basisobj);
fdnames   = getnames(harmfd);
x         = linspace(rangex(1), rangex(2), nx);
fdmat     = eval_fd(harmfd, x);
meanmat   = squeeze(eval_fd(pcastr.meanfd, x));
dimfd     = size(fdmat);
nharm     = dimfd(2);
casenames = getfdlabels(fdnames);
if harm == 0
    harm = (1:nharm);
end
if length(dimfd) == 2
    %  plotting for univariate functions
    for iharm = harm
        if matplt
            switch length(harm)
                case 1
                    subplot(1,1,1)
                case 2
                    subplot(2,1,iharm)
                case 3
                    subplot(3,1,iharm)
                case 4
                    subplot(2,2,iharm)
                otherwise
                    subplot(1,1,1)
            end
        else
            subplot(1,1,1)
        end    
        if expand == 0
            fac = sqrt(pcastr.values(iharm));
        else
            fac = expand;
        end
        vecharm    = fdmat(:,iharm);
        percentvar = round(100 * pcastr.varprop(iharm));
        meanplus   = meanmat+fac.*vecharm;
        meanminus  = meanmat-fac.*vecharm;
        plottop    = max([meanplus;meanminus]);
        plotbot    = min([meanplus;meanminus]);
        plot(x, meanmat,  '-')
        if colorwrd
            hold on
            plot(x, meanminus, 'r-')
            plot(x, meanplus,  'g-')
            hold off
        else
            text(x, meanplus,  '+')
            text(x, meanminus, '-')
        end
        xlabel(['\fontsize{12} ',fdnames{1}])
        ylabel(['\fontsize{12} Harmonic ',num2str(iharm)])
        axis([rangex(1), rangex(2), plotbot, plottop])
        if isempty(casenames)
            title(['PCA function ', num2str(iharm), ...
                ' (Percentage of variability ',  ...
                num2str(percentvar), ')'])
        else
            title(['PCA function ', casenames(iharm,:), ...
                ' (Percentage of variability ',  ...
                num2str(percentvar), ')'])
        end
        if matplt && length(harm) <= 4
            if length(harm) > 2 && length(harm) <= 4
                axis('square')
            end
        else
            disp('Press any key to continue')
            pause;
        end
    end
else
    %  plotting for multivariate functions
    if cycle && dimfd(3) == 2
        %  cycle plotting for  bivariate functions
        for iharm = harm
            if matplt
                switch length(harm)
                    case 1
                        subplot(1,1,1)
                    case 2
                        subplot(2,1,iharm)
                    case 3
                        subplot(3,1,iharm)
                    case 4
                        subplot(2,2,iharm)
                    otherwise
                        subplot(1,1,1)
                end
            else
                subplot(1,1,1)
            end
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
            xlabel(['\fontsize{12} ',fdnames{1}])
            ylabel(['\fontsize{12} Harmonic ',num2str(iharm)])
            axis([plotbot, plottop, plotbot, plottop])
            if isempty(casenames)
                title(['PCA function ', num2str(iharm), ...
                    ' (Percentage of variability ',  ...
                    num2str(percentvar), ')'])
            else
                title(['PCA function ', casenames(iharm,:), ...
                    ' (Percentage of variability ',  ...
                    num2str(percentvar), ')'])
            end
            if matplt && length(harm) <= 4
                if length(harm) > 2 && length(harm) <= 4
                    axis('square')
                end
            else
                disp('Press any key to continue')
                pause;
            end
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
            funnames = fdnames{3};
            for jvar = 1:nvar
                if iscell(funnames)
                    varnames = funnames{2};
                    varnamej = varnames(jvar,:);
                else
                    varnamej = ['Variable ',num2str(jvar)];
                end
                percentvar = round(100 * pcastr.varprop(iharm));
                vecharm    = fdmat(:,iharm,jvar);
                meanplus   = meanmat(:,jvar)+fac.*vecharm;
                meanminus  = meanmat(:,jvar)-fac.*vecharm;
                plottop    = max([meanplus;meanminus]);
                plotbot    = min([meanplus;meanminus]);
                subplot(nvar,1,jvar)
                plot(x, meanmat(:,jvar), 'b-', ...
                     x, meanplus,        'g-', ...
                     x, meanminus,       'r-')
                if jvar == nvar
                    xlabel(['\fontsize{12} ',fdnames{1}]);
                end
                ylabel(['\fontsize{12} ',varnamej])
                axis([rangex(1), rangex(2), plotbot, plottop])
                if jvar == 1
                    if isempty(casenames)
                        title(['\fontsize{13} Harmonic ',       ...
                               num2str(iharm),                  ...
                               ' (Percentage of variability ',  ...
                               num2str(percentvar), ')'])
                    else
                        title(['\fontsize{13} Harmonic ',       ...
                               casenames(iharm,:),              ...
                               ' (Percentage of variability ',  ...
                               num2str(percentvar), ')'])
                    end
                end
            end
            disp('Press any key to continue')
            pause;
        end        
    end
end

