function plot(basisobj, nx)
%  Plot a basis object.

%  last modified 15 June 2010

typex   = getbasistype(basisobj);
nbasisx = getnbasis(basisobj);

if ~strcmp(typex, 'FEM')
    
    %  set up fine mesh of values
    
    if nargin < 2, nx = max([10*nbasisx+1, 201]);  end
    
    %  evaluate basis at a fine mesh of values
    
    rangex   = getbasisrange(basisobj);
    x        = linspace(rangex(1),rangex(2),nx)';
    basismat = full(eval_basis(x, basisobj));
    
    %  plot the basis values
    
    phdl = plot (x, basismat, '-');
    set(phdl, 'LineWidth', 1);
    
    %  if the basis is of spline type, plot the knots
    
    if strcmp(typex, 'bspline')
        knots = getbasispar(basisobj);
        hold on
        for k=1:length(knots)
            lhdl = plot([knots(k), knots(k)], [0,1]);
            set(lhdl, 'LineWidth', 1, 'LineStyle', ':', 'color', 'r');
        end
        hold off
    end
    
    %  set plotting range
    
    if strcmp(typex, 'bspline')
        minval = 0;
        maxval = 1;
    else
        minval = min(min(basismat));
        maxval = max(max(basismat));
    end
    if minval == maxval
        if abs(minval) < 1e-1
            minval = minval - 0.05;
            maxval = maxval + 0.05;
        else
            minval = minval - 0.05*minval;
            maxval = maxval + 0.05*minval;
        end
    end
    xlabel('\fontsize{13} t')
    ylabel('\fontsize{13} \phi(t)')
    titstr = ['\fontsize{16} ', typex, ' basis', ...
        ',  no. basis fns = ', ...
        num2str(nbasisx)];
    if strcmp(typex, 'bspline')
        norderx = nbasisx - length(knots);
        titstr = [titstr, ',  order = ', num2str(norderx)];
    end
    title(titstr);
    axis([rangex(1), rangex(2), minval, maxval])
    
else
    params = getbasispar(basisobj);
    if isempty(params.dl)
        phdl=triplot(params.t,params.p(:,1),params.p(:,2));
        set(phdl, 'LineWidth', 2)
        hold on
        phdl=plot(params.nodes(:,1),params.nodes(:,2),'o');
        set(phdl, 'LineWidth', 2)
        phdl=plot(params.p(:,1),params.p(:,2),'.');
        set(phdl, 'LineWidth', 2)
        hold off
        xlabel('\fontsize{13} X')
        ylabel('\fontsize{13} Y')
        title('\fontsize{13} nodes at open o, vertices at filled o')
%         pdemesh(params.p', params.e', params.t')  
    else
        figure(1)
        pdegplot(dl) %  plot the decomposed geometry
        figure(2)
        pdemesh(params.p', params.e', params.t')  %  plot the mesh
    end
end
