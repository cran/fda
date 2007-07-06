function plotfit_fd(y, argvals, fd, casenames, varnames, residual, sortwrd, ...
                    rng, index, nfine)
%PLOTFIT plots discrete data along with a functional data object for fitting the
%  data.  It is designed to be used after DATA2FD, SMOOTH.FD or SMOOTH.BASIS to
%  check the fit of the data offered by the FD object.
%  Arguments:
%  Y         ... the data used to generate the fit
%  ARGVALS   ... discrete argument values associated with data
%  FD        ... a functional data object for fitting the data
%  CASENAMES ... a matrix of names for replicates
%  VARNAMES  ... a matrix of names for functions
%  RESIDUAL  ... if T, the residuals are plotted instead of the data plus curve
%  SORTWRD   ... sort plots by mean square error
%  RNG       ... a range of argument values to be plotted
%  INDEX     ... an index for plotting subsets of the curves (either sorted or not)
%  NFINE     ... number of points to use for plotting curves

%  Last modified 21 July 2006

if nargin < 10, nfine    = max([201,5*length(argvals)]); end
if nargin <  7, sortwrd  = 0;   end
if nargin <  6, residual = 0;   end

if size(argvals,1) == 1; argvals = argvals';  end

basis    = getbasis(fd);
rangeval = getbasisrange(basis);
if nargin < 8, rng = rangeval; end

coef  = getcoef(fd);
coefd = size(coef);
ndim  = length(coefd);

n  = size(y,1);
nrep = coefd(2);   
casenum = 1:nrep;
if ndim < 3, nvar = 1;  else nvar = coefd(3);  end
if nargin < 9, index = 1:nrep; end
y = reshape(y, n, nrep, nvar);
if nargin < 4, casenames = []; end

fdnames = getnames(fd);
argname = fdnames{1};
if nargin < 5 
    if nvar > 1
        varnames = []; 
    else
        varnames = fdnames{3}; 
    end
end

%  compute fitted values for evalargs and fine mesh of values

yhat   = reshape(eval_fd(argvals, fd),[n,nrep,nvar]);
res    = y - yhat;
MSE    = squeeze(mean(res.^2))';
MSE    = reshape(MSE,nvar,nrep);
MSEsum = sum(MSE);

%  compute fitted values for fine mesh of values

xfine = linspace(rng(1), rng(2), nfine)';
yfine = reshape(eval_fd(xfine, fd),[nfine,nrep,nvar]);

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
MSE = MSE  (:,index);
if ~isempty(casenames), casenames = casenames(index,:); end
casenum = casenum(index);
nrep  = length(index);

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
	        plot(argvals, res(:,i,j), '.', [rng(1),rng(2)], [0,0], ':')
            axis([rng(1),rng(2),ylimit(1),ylimit(2)])
	        xlabel(argname)
            if isempty(varnames)
                ylabel(['Residual for function ',num2str(j)])
            else
                ylabel(['Residual for ',varnames(j,:)])
            end
            if isempty(casenames)
	            title(['Case ',num2str(casenum(i)), ...
	                   '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
            else
	            title(['Case ',casenames(i,:), ...
	                   '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
            end
            if nrep > 1 || nvar > 1, pause;  end
       end
   end
 else
	%  plot the data and fit
	ylimit = [min(min(min(min(y))),min(min(min(yfine)))), ...
           max(max(max(max(y))),max(max(max(yfine))))];
    for i = 1:nrep 
        for j = 1:nvar
	        plot(argvals, y(:,i,j), '.', xfine, yfine(:,i,j), '-')
            axis([rng(1),rng(2),ylimit(1),ylimit(2)])
	        xlabel(argname)
            if isempty(varnames)
                ylabel(['Function ',num2str(j),' ',fdnames{3}])
            else
                ylabel(varnames(j,:))
            end
            if isempty(casenames)
 	            title(['Case ',num2str(casenum(i)), ...
	                       '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
            else
	            title(['Case ',casenames(i,:), ...
	                   '  RMS residual = ',num2str(sqrt(MSE(j,i)))])
            end
            if nrep > 1 || nvar > 1, pause;  end
        end
    end
end
