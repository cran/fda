function [pval,qval,Fobs,Fnull,Fvals,Fnullvals,pvals_pts,qvals_pts,fregresscell,argvals] = ...
    Fperm_fd(yfdPar, xfdlist, betalist,wt,nperm,argvals,q,plotres)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fperm_fd
%
% Computes a permutation F-test for the significance of a functional linear
% model, based on the maximal F-statistic for functional resposes, or an
% F-statistic for scalar responses. 
%
% Arguments:
%
%   - yfdPar, xfdlist, betalist, wt; arguments to fRegress for the
%   functional linear model. 
%
%   - nperm: number of permutations to use, defaults to 200
%
%   - argvals: time points to evaluate an F-statistic, if the response is
%   functional. Defaults to 101 equally spaced points on the range of the
%   response. 
% 
%   - q: alpha level of the test, defaults to 0.05.
%
%   - plotres: should a graphical display be given?
%
% Returns:
%   - pval: the observed p-value for the test
%   
%   - qval: the permutation critical value
%
%   - Fobs: the observed test statistic: maximal pointwise F-statistic
%   
%   - Fnull: the computed permutation distribution given as a vector of
%   samples
%
%   - Fvals: the pointwise t-statistics
%
%   - Fnullvals: a matrix of pointwise null distribution values
%
%   - qvals_pts: the pointwise critical values for the permutation
%   distribution
%
%   - pvals_pts: pointwise p-values for the permutation test
%
%   - fregresscell: the functional linear model as returned from fRegress. 
%
%   - argvals: the time points at which the t-statistics were evaluated.

% Last modified 15 May, 2009

   if nargin<8, plotres = 1; end
   if nargin < 7, q = 0.05; end
   if nargin < 6, argvals = []; end
   if nargin < 5, nperm = 200; end
    
  Fnull = zeros(nperm,1);
  Fnullvals = [];

  q = 1-q;

  tic;
  fregresscell = fRegress(yfdPar, xfdlist, betalist);
  elapsed.time = toc;

  if( elapsed.time > 30/nperm )
    disp(strcat('Estimated Computing time =  ',num2str(round(nperm*elapsed.time)),' seconds.'))
  end

    yhat = fregresscell{5};

%%  if isa_fdPar(fregresscell.yfdPar), yfdPar = getfd(fregresscell.fdPar); 
%%  else yfdPar = fregresscell.yfdPar

  [Fvals,argvals] = Fstat_fd(yfdPar,yhat,argvals);

  Fobs = max(Fvals);

  if isnumeric(yfdPar), n = length(yfdPar);
  else
      tempc = getcoef(yfdPar);
      n = size(tempc,2); 
  end

  for i = 1:nperm

    tyfdPar = yfdPar(randperm(n));

    tfregresscell = fRegress(tyfdPar, xfdlist, betalist);
    yhat = tfregresscell{5};
    
    Fnullvals = [Fnullvals,Fstat_fd(yfdPar,yhat,argvals)];

    Fnull(i) = max(Fnullvals(:,i));
  end


    pval = mean( Fobs < Fnull );
    qval = quantile(Fnull,q);

    pvals_pts = mean(repmat(Fvals,1,nperm)<Fnullvals,2);
    qvals_pts = quantile(Fnullvals,q,2);


    
    if plotres 
        if isa_fd(yfdPar),
            fdnames = getnames(yfdPar);
            
            plot(argvals,Fvals,'r','linewidth',2);
			xlabel(fdnames(3));
            ylabel('F-statistic');
            title('Permutation F-Test');
            hold on
            plot(argvals,qvals_pts,'b-.','linewidth',2);
            plot([min(argvals) max(argvals)],qval*ones(1,2),'b-','linewidth',2);
            legend('Observed Statistic',strcat('pointwise ',num2str(1-q),' critical value'),strcat('maximum ',num2str(1-q),' critical value'));
        else
            N = hist(Fnull);
            xlabel('F-value')
			title('Permutation F-Test')
            hold on
            plot(Fobs*ones(1,2),[0 max(N)],'r-')
            plot(qval*ones(1,2),[0 max(N)],'b-')
            legend('Observed Statistic',strcat('Permutation ',1-q,' critical value'))

        end
    end
