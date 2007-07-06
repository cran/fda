function plot(fdobj, Lfdobj, matplt, href, nx)
%  PLOT   Plot a functional data object.
%  Arguments:
%  FDOBJ   ... A functional data object to be plotted.
%  LFDOBJ  ... A linear differential operator object or a positive
%              integer specifying a derivative.  This operator is
%              applied to the functions before plotting.
%  MATPLOT ... If MATPLT is nonzero, all curves are plotted in a 
%              single plot.
%              Otherwise, each curve is plotted separately, and the
%              next curve is plotted when the mouse is clicked.
%  HREF    ... If HREF is nonzero, a horizontal dotted line through 0 
%              is plotted.
%  NX      ... Number of equally spaced argument values to use for 
%              plotting.  Default is 101.  

%  Last modified 20 July 2006

%  set default arguments

if nargin < 5, nx = 101;             end
if nargin < 4, href = 1;             end
if nargin < 3, matplt = 1;           end
if nargin < 2, Lfdobj = int2Lfd(0);  end

%  check arguments

if ~isa_fd(fdobj)
    error ('Argument fdobj is not a functional data object.');
end

Lfdobj = int2Lfd(Lfdobj);
if ~isa_Lfd(Lfdobj)
    error ('Argument Lfdobj is not a linear differential operator object.');
end

% set up dimensions of problem

coef    = getcoef(fdobj);
coefd   = size(coef);
ndim    = length(coefd);
ncurve  = coefd(2);

if ndim > 2
    nvar = coefd(3);
else
    nvar = 1;
end

basisobj = getbasis(fdobj);
rangex   = getbasisrange(basisobj);
x        = linspace(rangex(1),rangex(2),nx)';
fdmat    = eval_fd(x, fdobj, Lfdobj);

fdnames  = getnames(fdobj);

switch ndim
    case 1
        frng(1) = min(fdmat);
        frng(2) = max(fdmat);
    case 2
        frng(1) = min(min(fdmat));
        frng(2) = max(max(fdmat));
    case 3
        frng(1) = min(min(min(fdmat)));
        frng(2) = max(max(max(fdmat)));
    otherwise
        frng = [1 1];
end

if frng(1) == frng(2)
    if abs(frng(1)) < 1e-1
        frng(1) = frng(1) - 0.05;
        frng(2) = frng(2) + 0.05;
    else
        frng(1) = frng(1) - 0.05*abs(frng(1));
        frng(2) = frng(2) + 0.05*abs(frng(1));
    end
end

if ndim < 2
    if href && (frng(1) <= 0 && frng(2) >= 0)
        plot(x,fdmat, '-', x, zeros(nx), ':');
    else
        plot (x,fdmat, '-');
    end
    xlabel(fdnames{1});
    ylabel(fdnames{3});
    axis([x(1), x(nx), frng(1), frng(2)])
end

if ndim == 2
    if matplt
        if href && (frng(1) <= 0 && frng(2) >= 0)
            plot (x, fdmat, '-', x, zeros(nx), ':');
        else
            plot (x, fdmat, '-');
        end
        xlabel(fdnames{1});
        ylabel(fdnames{3});
        if frng(2) > frng(1), axis([x(1), x(nx), frng(1), frng(2)]); end
    else
        for icurve = 1:ncurve
            plot (x, fdmat(:,icurve), '-');
            title(['Curve ', num2str(icurve)]);
            if href && (frng(1) <= 0 && frng(2) >= 0)
                plot(x, zeros(nx), ':');
            end
            xlabel(fdnames{1});
            ylabel(fdnames{3});
            if frng(2) > frng(1), axis([x(1), x(nx), frng(1), frng(2)]); end
            ginput(1);
        end
    end
end

if ndim == 3
    if matplt
        for ivar = 1:nvar
            subplot(nvar,1,ivar);
            temp = squeeze(fdmat(:,:,ivar));
            if href && (frng(1) <= 0 && frng(2) >= 0)
                plot (x, temp, '-', x, zeros(nx), ':');
            else
                plot (x, temp, '-');
            end
            xlabel(fdnames{1});
            ylabel(fdnames{3});
            axis([x(1), x(nx), frng(1), frng(2)])
            title(['Variable ', num2str(ivar)]);
        end
    else
        for icurve = 1:ncurve
            for ivar = 1:nvar
                subplot(nvar,1,ivar);
                temp = squeeze(fdmat(:,icurve,ivar));
                if href && (frng(1) <= 0 && frng(2) >= 0)
                    plot (x, temp, '-', x, zeros(nx), '.');
                else
                    plot (x, temp, '-');
                end
                xlabel(fdnames{1});
                ylabel(fdnames{3});
                axis([x(1), x(nx), frng(1), frng(2)])
                title(['Variable ', num2str(ivar), ' Curve ', num2str(icurve)]);
            end
            ginput(1);
        end
    end
end

