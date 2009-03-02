function [pval,qval,Tobs,Tnull,Tvals,Tnullvals,qvals_pts,pvals_pts,argvals] = tperm_fd(x1fd,x2fd,nperm,q,argvals,plotres) 
%%%%%%%%%%%%%%%%%%%%%%%
% tperm_fd 
%
% Permutation t-tests between groups of functional data objects.
% Specifically, the function reports the permutation test for the maximal
% pointwise t-statistic. It also reports the permutation critical values
% for the t-statistic at each time point. 
%
% Arguments
%    - x1fd, x2fd  two functional data objects, replications within the fd
%    objects are the data of interest and the test is for differences in
%    the means of x1fd and x2fd.
%   
%   - nperm:  number of permutations to use, defaults to 200
%
%   - q: alpha-level of the test, defaults to 0.05
%
%   - argvals: number of points at which to compute the t-statistic.
%   Defaults to 100 equally spaced along the range of the functional data
%   objects. 
%
%   - plotres: should a graphical desplay be plotted?
%
% Returns
%   - pval: the observed p-value for the test
%   
%   - qval: the permutation critical value
%
%   - Tobs: the observed test statistic: maximal pointwise t-statistic
%   
%   - Tnull: the computed permutation distribution given as a vector of
%   samples
%
%   - Tvals: the pointwise t-statistics
%
%   - Tnullvals: a matrix of pointwise null distribution values
%
%   - qvals_pts: the pointwise critical values for the permutation
%   distribution
%
%   - pvals_pts: pointwise p-values for the permutation test
%
%   - argvals: the time points at which the t-statistics were evaluated.

% Last modified 15 May, 2009

    if ~isa_fd(x1fd) || ~isa_fd(x2fd)                                     
        error('x1fd and x2fd must both be functional data objects');         
    end                                                                      
                                                                           
    rangeobs = getbasisrange(getbasis(x1fd));
    rangehat = getbasisrange(getbasis(x2fd));


    if any(rangeobs ~= rangehat)
        error('x1fd and x2fd do not have the same range.')
    end

    if nargin < 6, plotres = 1; end
    if nargin < 5, argvals = linspace(rangeobs(1),rangeobs(2),101); end
    if nargin < 4, q = 0.05; end
    if nargin < 3, nperm = 200; end 
    

    q = 1-q;

    x1mat = eval_fd(argvals,x1fd);
    x2mat = eval_fd(argvals,x2fd);

    n1 = size(x1mat,2);
    n2 = size(x2mat,2);

    Xmat = [x1mat x2mat];

    Tnull = zeros(nperm,1);

    Tnullvals = zeros(length(argvals),nperm);

    for i = 1:nperm
        tXmat = Xmat(:,randperm(n1+n2));

        tmean1 = mean(tXmat(:,1:n1),2);
        tmean2 = mean(tXmat(:,n1+(1:n2)),2);

        tvar1 = var(tXmat(:,1:n1),[],2)/n1;
        tvar2 = var(tXmat(:,n1+(1:n2)),[],2)/n2;

        Tnullvals(:,i) = abs(tmean1-tmean2)./sqrt(tvar1+tvar2);
        Tnull(i) = max(Tnullvals(:,i));
    end

    mean1 = mean(Xmat(:,1:n1),2);
    mean2 = mean(Xmat(:,n1+(1:n2)),2);

    var1 = var(Xmat(:,1:n1),[],2)/n1;
    var2 = var(Xmat(:,n1+(1:n2)),[],2)/n2;

    Tvals = abs(mean1-mean2)./sqrt(var1+var2);
    Tobs = max(Tvals);

    pval = mean( Tobs < Tnull );
    qval = quantile(Tnull,q);
    
    pvals_pts = mean(repmat(Tvals,1,nperm)<Tnullvals,2);
    qvals_pts = quantile(Tnullvals,q,2);

    if plotres

        fdnames = getnames(x1fd);

        plot(argvals,Tvals,'r','linewidth',2)
        ylabel('t-statistic')
        xlabel(fdnames(3))
        hold on
        plot(argvals,qvals_pts,'r-.','linewidth',2);
        plot(rangeobs,qval*ones(1,2),'b-','linewidth',2);
        hold off
        
	
        legend('Observed Statistic',strcat('pointwise ',1-q,' critical value'),...
			    strcat('maximum ',1-q,' critical value'));

    end