function plotbeta(betaestcell, betastderrcell, argvals)
%  PLOTBETA plots a functional parameter along with confidence
%  limits
%  Arguments:
%  BETAESTCELL    ... A cell object containing one or more functional
%                     parameter objects or functional data objects.
%  BETASTDERRCELL ... A cell object containing functional data objects
%                     for the standard error of the objects in
%                     BETAESTCELL.

%  Last modified 20 July 2006

%  check BETAESTCELL

if strcmp(class(betaestcell), 'fdPar') || ...
   strcmp(class(betaestcell), 'fd')
    temp{1} = betaestcell;
    betaestcell = temp;
end

if ~iscell(betaestcell)
    error('BETAESTCELL is not a cell, fd, or fdpar object.');
end

%  check BETASTDERRCELL

if nargin > 1

    if strcmp(class(betastderrcell), 'fd')
        temp{1} = betastderrcell;
        betastderrcell = temp;
    end

    if ~iscell(betastderrcell)
        error('BETASTDERRCELL is not a cell, or fd object.');
    end

end

%  get range

if     isa_fdPar(betaestcell{1})
    rangeval = getbasisrange(getbasis(getfd(betaestcell{1})));
elseif isa_fd(betaestcell{1})
    rangeval = getbasisrange(getbasis(betaestcell{1}));
else
    error(['A cell does not contain either a functional parameter ', ...
           'or a functional data object.']);
end

if nargin < 3
    argvals = linspace(rangeval(1),rangeval(2),51)';
end
n = length(argvals);
p = length(betaestcell);
for j=1:p
    if     isa_fdPar(betaestcell{j})
        betavec = eval_fd(argvals, getfd(betaestcell{j}));
    elseif isa_fd(betaestcell{j})
        betavec = eval_fd(argvals, betaestcell{j});
    else
        error(['A cell of BETAESTCELL ', ...
               'does not contain either a functional parameter ', ...
               'or a functional data object.']);
    end
    plot(argvals, betavec, '-')
    line([rangeval(1),rangeval(2)],[0,0],'LineStyle',':','Color','r')
    if nargin > 1
        betastderr = eval_fd(argvals, betastderrcell{j});
        betavecp = betavec + 2.*betastderr;
        betavecm = betavec - 2.*betastderr;
        fill([argvals;   argvals(n:-1:1)], ...
             [betavecp;   betavecm(n:-1:1)], 'w', 'LineStyle', '--')
        line(argvals, betavec,'Color','b','LineWidth',2)
        line([rangeval(1),rangeval(2)],[0,0],'LineStyle',':','Color','r')
        for i=1:n
            line([argvals(i),argvals(i)],[betavecm(i),betavecp(i)])
        end
    end
    title(['\fontsize{16} Regression function ',num2str(j)])
    v = axis;
    axis([rangeval(1),rangeval(2),v(3),v(4)])
    if p > 1, pause; end
end
