function plot_FEM(fdobj, X, Y, nderivs, nfine)
% PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by 
% vectors X and Y;
%
% Last modified on 15 May 2011

if nargin < 5,                      nfine = 101;      end
if nargin < 4 || isempty(nderivs),  nderivs = [0,0];  end
if nargin < 3,                      Y = [];           end
if nargin < 2,                      X = [];           end

if ~isa_fd(fdobj)
   error('FDOBJ is not an FD object')
end

%  check derivatives

if length(nderivs) ~= 2
    error('NDERIVS not of length 2.');
end
if sum(nderivs) > 2
    error('Maximum derivative order is greater than two.');
end

coefmat = full(getcoef(fdobj));
nsurf = size(coefmat,2);

basisobj = getbasis(fdobj);

params = getbasispar(basisobj);

p         = params.p;
t         = params.t;
t         = t(:,1:3);

if isempty(X)
    xmin = min(p(:,1));
    xmax = max(p(:,1));
    nx   = nfine;
    X    = linspace(xmin, xmax, nx)';
else
    nx   = length(X);
end

if isempty(Y)
    ymin = min(p(:,2));
    ymax = max(p(:,2));
    ny   = nfine;
    Y    = linspace(ymin, ymax, ny)';    
else
    ny   = length(Y);
end

Xmat = X*ones(1,ny);
Ymat = ones(nx,1)*Y';
Xvec = Xmat(:);
Yvec = Ymat(:);


evalmat = zeros(nx*ny,nsurf);
for iopr = 1:size(nderivs,1)
    evalmat = evalmat + eval_FEM_fd(Xvec, Yvec, fdobj, nderivs(iopr,:));
end

for isurf=1:nsurf
    evalmati = reshape(evalmat(:,isurf),nx,ny)';
    surf(X,Y,evalmati);
    xlabel('\fontsize{13} X')
    ylabel('\fontsize{13} Y');
    if nsurf > 1
        pause
    end
end


