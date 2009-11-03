function Fpermstr = Fperm_fd(yfdPar, xfdcell, betacell, wt, ...
                             nperm, argvals, q, plotres)
%  FPERM_FD does a permutation test for a functional parameter effect.
%  Arguments;
%  The first four are same as the corresponding arguments for function
%  fRegress; yfdpar, xfdCell, betacell, wt
%  NPERM   ... number of permutations to be used
%  ARGVALS ... argument values for evaluating functional responses,
%  Q       ... tail probability of quantile to compare
%  PLOTRES ... If nonzero, results are plotted
%  Results;
%  FPERMSTR;
%  A struct object with slots;
%

%  Last modified 22 September 2009 by Jim Ramsay

%  Set default arguments

if nargin < 8, plotres = 1;  end
if nargin < 7, q = 0.05;     end
if nargin < 6, argvals = []; end
if nargin < 5, nperm = 200;  end
if nargin < 4, wt = [];      end

%

q = 1-q;

if isa_fd(xfdcell{1})
    N = size(getcoef(xfdcell{1}),2);
else
    N = length(xfdcell{1});
end

tic;
fRegressStr = fRegress(yfdPar, xfdcell, betacell, wt);
elapsed_time = toc;

if elapsed_time > 30/nperm
    disp(['Estimated Computing time = ', ...
        num2str(nperm*elapsed_time),' seconds.'])
end

yhatfd = fRegressStr.yhat;

[Fvals, argvals] = Fstat_fd(yfdPar,yhatfd,argvals);
nargs = length(argvals);

Fobs = max(Fvals);

Fnull     = zeros(nperm,1);

Fnullvals = [];

for i = 1:nperm
    yfdPari      = yfdPar(randperm(N));
    fRegressStri = fRegress(yfdPari, xfdcell, betacell, wt);
    yhatfdi      = fRegressStri.yhat;
    Fi           = Fstat_fd(yfdPar,yhatfdi,argvals);
    Fnullvals    = [Fnullvals,Fi];
    Fnull(i) = max(Fnullvals(:,i));
end

pval = mean(Fobs < Fnull);
qval = quantile(Fnull,q);

pvals_pts = zeros(nargs,1);
qvals_pts = pvals_pts;
for i=1:nargs
    pvals_pts(i) = mean(Fvals(i) < Fnullvals(:,i));
    qvals_pts(i) = quantile(Fnullvals(:,i),q);
end

if plotres
    if isa_fd(yfdPar)
        plot(argvals, Fvals, 'b-', ...
            argvals, qvals_pts, 'b--', ...
            [min(argvals),max(argvals)], [qval, qval], 'b:')
        xlabel('\fontsize{13} argument values')
        ylabel('\fontsize{13} F-statistic')
        title('\fontsize{13} Permutation F-Test')
        legend('Observed Statistic', ...
            ['pointwise ',num2str(1-q),' critical value'], ...
            ['maximum ',  num2str(1-q),' critical value'])
    else
        cnts  = hist(Fnull);
        xvals = [Fnull;Fobs];
        xmax  = ceil(max(xvals));
        xmin  = min(xvals);
        Nmax  = max(cnts);
        hist(Fnull)
        hold on
        phdl=plot([Fobs, Fobs], [0,Nmax], 'b-', ...
                  [qval, qval], [0,Nmax], 'b--');
        set(phdl, 'LineWidth', 2)
        hold off
        axis([xmin,xmax,0, max(N)])
        xlabel('\fontsize{13} F-value')
        title('\fontsize{13} Permutation F-Test')
        legend('Frequency', 'Observed Statistic', ...
            ['Permutation ',num2str(1-q),' critical value'])        
    end
end

Fpermstr.pval        = pval;
Fpermstr.qval        = qval;
Fpermstr.Fobs        = Fobs;
Fpermstr.Fnull       = Fnull;
Fpermstr.Fvals       = Fvals;
Fpermstr.Fnullvals   = Fnullvals;
Fpermstr.pvals_pts   = pvals_pts;
Fpermstr.qvals_pts   = qvals_pts;
Fpermstr.fRegressStr = fRegressStr;
Fpermstr.argvals     = argvals;

