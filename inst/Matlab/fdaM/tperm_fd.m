function tpermStr = tperm_fd(x1fd, x2fd, nperm, q, argvals, plotres)
%  TPERM computes a permutation t-test of the difference between
%  the functional means of two groups of data.
%
%  Arguments:
%  X1FD    ... First  sample of functional observations
%  X2FD    ... Second sample of functional observations
%  NPERM   ... Number of permutations ... default 200
%  Q       ... Tail probability for test ... default 0.05
%  ARGVALS ... Argument values at which to evaluate functions
%  PLOTRES ... If nonzero, plot results ... default 0
%  Results:

%  Last modified 10 October 2009 by Jim Ramsay

%  check first two arguments

if nargin < 2
    error('Arguments X1FD and X2FD not supplied.');
end

if ~isa_fd(x1fd) || ~isa_fd(x2fd)
    error('x1fd and x2fd must both be functional data objects');
end

if length(size(getcoef(x1fd))) > 2 || length(size(getcoef(x2fd))) > 2
    error('Both of X1FD and X2FD are not univariate.');
end

%  Set default arguments

if nargin < 6,  plotres = 1;   end
if nargin < 5,  argvals = [];  end
if nargin < 4,  q = 0.05;      end
if nargin < 3,  nperm = 200;   end


range1 = getbasisrange(getbasis(x1fd));
range2 = getbasisrange(getbasis(x2fd));


if ~all(range1 == range2)
    error('x1fd and x2fd do not have the same range.');
end

if isempty(argvals)
    narg = 101;
    argvals = linspace(range1(1),range1(2),narg)';
else
    narg = length(argvals);
end

q = 1-q;

x1mat = eval_fd(argvals,x1fd);
x2mat = eval_fd(argvals,x2fd);

Xmat = [x1mat,x2mat];

n1 = size(x1mat,2);
n2 = size(x2mat,2);

Tnull = zeros(nperm,1);

Tnullvals = zeros(length(argvals),nperm);

for i = 1:nperm
    tXmat = Xmat(:,randperm(n1+n2));
    
    tmean1 = mean(tXmat(:,1:n1),     2);
    tmean2 = mean(tXmat(:,n1+(1:n2)),2);
    
    tvar1 = var(tXmat(:,1:n1),     0,2)/n1;
    tvar2 = var(tXmat(:,n1+(1:n2)),0,2)/n2;
    
    Tnullvals(:,i) = abs(tmean1-tmean2)./sqrt(tvar1+tvar2);
    Tnull(i) = max(Tnullvals(:,i));
end

mean1 = mean(Xmat(:,1:n1),     2);
mean2 = mean(Xmat(:,n1+(1:n2)),2);

var1 = var(Xmat(:,1:n1),     0,2)/n1;
var2 = var(Xmat(:,n1+(1:n2)),0,2)/n2;

Tvals = abs(mean1-mean2)./sqrt(var1+var2);
Tobs  = max(Tvals);

pval = mean( Tobs < Tnull );
qval = quantile(Tnull, q);

pvals_pts = zeros(narg,1);
qvals_pts = pvals_pts;
for i=1:narg
    pvals_pts(i) = mean(Tvals(i) < Tnullvals(i,:));
    qvals_pts(i) = quantile(Tnullvals(i,:), q);
end

if plotres
    
    ylim = [ min([Tvals;qvals_pts]),max([Tobs;qval])];
    
    phdl=plot(argvals, Tvals, 'b-', argvals, qvals_pts, 'b--', ...
        [min(argvals),max(argvals)], [qval,qval], 'r:');
    set(phdl, 'LineWidth', 2)
    xlabel('\fontsize{13} argument values')
    ylabel('\fontsize{13} t-statistic')
    axis([min(argvals),max(argvals),ylim])
    legend('\fontsize{13} Observed Statistic', ...
           ['\fontsize{13} pointwise ', 1-q, ' critical value'], ...
           ['\fontsize{13} maximum ',   1-q, ' critical value'], ...
           'location', 'NorthWest')  
end

tpermStr.pval         = pval;
tpermStr.qval         = qval;
tpermStr.Tobs         = Tobs;
tpermStr.Tnull        = Tnull;
tpermStr.Tvals        = Tvals;
tpermStr.Tnullvals    = Tnullvals;
tpermStr.pvals_pts    = pvals_pts;
tpermStr.qvals_pts    = qvals_pts;
tpermStr.argvals      = argvals;

