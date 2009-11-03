function plotfit_fd(y, argvals, fdobj, residual, sortwrd, rng, index)
%PLOTFIT plots discrete data along with a functional data object for 
%  fitting the data.  It is designed to be used after SMOOTH.BASIS to
%  check the fit of the data offered by the FD object.
%  Note:  As of March 2, 2009, the arguments CASENAMES, VARNAMES and
%  NFINE have been removed.
%  Arguments:
%  Y         ... the data used to generate the fit
%  ARGVALS   ... discrete argument values associated with data
%  FD        ... a functional data object for fitting the data
%  RESIDUAL  ... if nonzero, the residuals are plotted instead of 
%                the data plus curve
%  SORTWRD   ... sort plots by mean square error
%  RNG       ... a range of argument values to be plotted
%  INDEX     ... an index for plotting subsets of the curves 
%                (either sorted or not)
%  NFINE     ... number of points to use for plotting curves

%  Last modified 24 December 2011

if nargin < 5 || isempty(sortwrd),  sortwrd  = 0;   end
if nargin < 4 || isempty(residual), residual = 0;   end

if size(argvals,1) == 1; argvals = argvals';  end

basisobj = getbasis(fdobj);
rangeval = getbasisrange(basisobj);
nbasis   = getnbasis(basisobj);

if nargin < 8 || isempty(rng), rng = rangeval; end

coef  = getcoef(fdobj);
coefd = size(coef);
ndim  = length(coefd);

n    = size(y,1);
nrep = coefd(2);   
casenum = 1:nrep;
if ndim < 3, nvar = 1;  else nvar = coefd(3);  end

if nargin < 6 || isempty(index), index = 1:nrep;  end

y = reshape(y, n, nrep, nvar);

%  set up number of points at which to evaluate the curves

nfine = max([201, 10*nbasis+1]);

%  set up labels for arguments, cases and variables.

fdnames   = getnames(fdobj);
argname   = fdnames{1};
casenames = fdnames{2};
varnames  = fdnames{3};

%  compute fitted values for evalargs and fine mesh of values

yhat   = reshape(eval_fd(argvals, fdobj),[n,nrep,nvar]);
res    = y - yhat;
MSE    = squeeze(mean(res.^2))';
MSE    = reshape(MSE,nvar,nrep);
MSEsum = sum(MSE);

%  compute fitted values for fine mesh of values

xfine = linspace(rng(1), rng(2), nfine)';
yfine = reshape(eval_fd(xfine, fdobj),[nfine,nrep,nvar]);

%  sort cases by MSE if desired

if sortwrd && nrep > 1
    if nvar > 1
	    [MSEj,MSEind] = sort(MSEsum);
        MSE = MSE(:,MSEind);
    else
        [MSE,MSEind]  = sort(MSE);
    end
	y      = y(:,MSEind,:);
	yfine  = yfine(:,MSEind,:);
	res    = res  (:,MSEind,:);
    if ~isempty(casenames), casenames = casenames(MSEind,:); end
    casenum = casenum(MSEind);
end

%  set up fit and data as 3D arrays, selecting curves in INDEX

y     = y    (:,index,:);
res   = res  (:,index,:);
yfine = yfine(:,index,:);
MSE   = MSE  (:,index);
% if ~isempty(casenames), casenames = casenames(index,:); end
% casenum = casenum(index);
nrep    = length(index);

%  select values in ARGVALS, Y, and YHAT within RNG

argind  = argvals >= rng(1) & argvals <= rng(2);
argvals = argvals(argind);
y       = y   (argind,:,:);
res     = res (argind,:,:);

xfiind = xfine >= rng(1) & xfine <= rng(2);
xfine  = xfine(xfiind);
yfine  = yfine(xfiind,:,:);
	
%  plot the results

if residual
	%  plot the residuals
	ylimit = [min(min(min(res))), max(max(max(res)))];
    for i = 1:nrep 
        for j = 1:nvar
            subplot(nvar,1,j)
	        plot(argvals, res(:,i,j), '.', [rng(1),rng(2)], [0,0], ':')
            axis([rng(1),rng(2),ylimit(1),ylimit(2)])
	        xlabel(['\fontsize{12} ',argname])
            if iscell(varnames)
                ylabel(['\fontsize{12} ', varnames{2}(j,:), ' residual'])
            else
                ylabel(['\fontsize{12} ', varnames, ' residual'])
            end
            if nrep > 1
                if iscell(casenames)
                    title(['\fontsize{12} ', casenames{2}(j,:), ...
                        '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
                else
                    title(['\fontsize{12} Case ', num2str(i), ...
                        '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
                end
            else
                if iscell(casenames)
                    title(['\fontsize{12} RMS residual = ',num2str(sqrt(MSE(j,i)))])
                else
                    title(['\fontsize{12} RMS residual = ',num2str(sqrt(MSE(j,i)))])
                end
            end
        end
        pause
   end
 else
	%  plot the data and fit
	ylimit = [min(min(min(min(y))),min(min(min(yfine)))), ...
           max(max(max(max(y))),max(max(max(yfine))))];
    for i = 1:nrep 
        for j = 1:nvar
	        phdl = plot(argvals, y(:,i,j), '.', xfine, yfine(:,i,j), 'b-');
            set(phdl, 'LineWidth', 1.5)
            axis([rng(1),rng(2),ylimit(1),ylimit(2)])
	        xlabel(['\fontsize{12} ',argname])
            if iscell(varnames)
                ylabel(['\fontsize{12} ', varnames{2}(j,:)])
            else
                ylabel(['\fontsize{12} ', varnames])
            end
            if nrep > 1
                if iscell(casenames)
                    title(['\fontsize{12} ', casenames{2}(j,:), ...
                        '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
                else
                    title(['\fontsize{12} Case ', num2str(i), ...
                        '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
                end
            else
                if iscell(casenames)
                    title(['\fontsize{12} RMS residual = ',num2str(sqrt(MSE(j,i)))])
                else
                    title(['\fontsize{12} RMS residual = ',num2str(sqrt(MSE(j,i)))])
                end
            end
            if nrep > 1 || nvar > 1, pause;  end
        end
    end
end
