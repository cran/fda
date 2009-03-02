function yregfd = register_newfd(yfd, Wfd,type)

if nargin < 3, type = 'direct'; end

coef  = getcoef(Wfd);
shift = coef(1);

if all(coef == 0)
    if strcmp(type,'periodic'),
        if shift == 0
            yregfd = yfd;
        end
    else
        yregfd = yfd;
    end
else

    % Now evaluate on a fine grid:

    Wrange = getbasisrange(getbasis(Wfd));
    ybasis = getbasis(yfd);
    yrange = getbasisrange(ybasis);

    if strcmp(type,'periodic') && any(Wrange ~= yrange)
        error('Registration functions and functions to be registered must have the same range')
    end

    neval  = max(10*getnbasis(ybasis) + 1,101);
    tfine  = linspace(yrange(1),yrange(2),neval);

    % Shift registration is easy

    if strcmp(type,'periodic'),
        yfine = eval_fd(tfine,yfd);
        yfine = shifty(tfine,yfine,shift);
        ycoef  = project_basis(yfine, tfine, ybasis, 1);
        yregfd = fd(ycoef, ybasis);
    end

    % On the other hand, if we have warping functions:

    if  strcmp(type,'direct'), xfine = eval_fd(tfine,Wfd); end

    if strcmp(type,'monotone'),
        xfine = eval_monfd(tfine,Wfd);
        xfine = xfine*diag(1./xfine(neval,:))*(Wrange(2)-Wrange(1))+Wrange(1);
    end
    xfine = xfine.*( xfine>yrange(1) & xfine < yrange(2)) + yrange(2).*(xfine>=yrange(2)) + yrange(1).*(xfine<=yrange(1));
    yfine = eval_fd(tfine,yfd);

    xdim = size(xfine);
    ydim = size(yfine);

    % Check that we have the right dimensions, we can register multiple y dimensions
    % to one warping function, but we must have as many warping functions as there
    % are y replicates

    if xdim(2) ~= ydim(2),
        error('There must be as many warping function replicates as y replicates')
    end
    if length(ydim) == 3 && length(xdim)==2, xfine = array(xfine,ydim); end


    % Now do the registration

    ycoef = 0*getcoef(yfd);
    cdim = size(ycoef);

    if length(cdim)==1
        yfine = eval_fd(xfine,yfd);
        ycoef = project_basis(yfine,tfine,ybasis,1);
    end

    if length(cdim)==2
        for i = 1:cdim(2)
            yfine = eval_fd(xfine(:,i),yfd(i));
            ycoef(:,i) = project_basis(yfine,tfine,ybasis);
        end
    end

    if length(cdim)==3
        for j = 1:cdim(3)
            for i = 1:cdim(2)
                yfine = eval_fd(xfine(:,i,j),yfd(i,j));
                ycoef(:,i,j) = project_basis(yfine,tfine,ybasis);
            end
        end
    end

    yregfd = fd(ycoef, ybasis);

end