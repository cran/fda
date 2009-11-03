addpath ('c:\Matlab7\fdaM')

load herm1.dat
load herm2.dat
load herm3.dat
load herm4.dat
load herm5.dat
load herm6.dat
load herm7.dat

%  ---------------------------------------------------------
%               plot results for the 6th baby
%  ---------------------------------------------------------

x = herm6(:,1);
y = herm6(:,2);

plot(x, y, 'o')

n  = length(x);
wt = ones(n,1);
zmat = wt;

rng = [1,n];

nbasis = n + 2;
basis  = create_bspline_basis([1,n], nbasis);
Wfd0   = fd(zeros(nbasis,1), basis);
WfdPar = fdPar(Wfd0, 2, 1e-4);
[Wfd, beta] = smooth_monotone(x, y, WfdPar, zmat, wt);

xfine = linspace(1,n,151)';
yhat  = beta(1) + beta(2).*eval_mon(xfine,Wfd);
Dyhat =           beta(2).*eval_mon(xfine, Wfd, 1);

ahdl = axes('Box', 'on', 'FontSize', 16);
set(ahdl, 'FontSize', 16)
set(ahdl, 'Xlim', [0,40]);
set(ahdl, 'Xtick', [0:10:40]);
set(ahdl, 'Ylim', [112,133]);
set(ahdl, 'Ytick', 115:5:130);
xlabel('Day', 'FontSize', 19);
ylabel('Tibia length (mm)', 'FontSize', 19);
lhdl = line(xfine, yhat);
set(lhdl, 'LineWidth', 2)
lhdl = line(x, y);
set(lhdl, 'LineWidth', 1, 'LineStyle', 'none', 'Marker', 'o')

print -dps2 'c:/MyFiles/fdabook/revision/figs.dir/growonebabyH.ps'

ahdl = axes('Box', 'on', 'FontSize', 16);
set(ahdl, 'FontSize', 16)
set(ahdl, 'Xlim', [0,40]);
set(ahdl, 'Xtick', [0:10:40]);
set(ahdl, 'Ylim', [0,2.3]);
set(ahdl, 'Ytick', 0:0.5:2);
xlabel('Day', 'FontSize', 19);
ylabel('D Tibia length (mm/day)', 'FontSize', 19);
lhdl = line(xfine, Dyhat);
set(lhdl, 'LineWidth', 2)

print -dps2 'c:/MyFiles/fdabook/revision/figs.dir/growonebabyV.ps'




