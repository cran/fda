function evalarray = eval(bifd, sevalarg, tevalarg, sLfd, tLfd)
%EVAL  evaluates a bi-functional data object BIFD
%  at argument values in arrays SEVALARG and TEVALARG.
%  SLfd and TLfd are either integers giving the order of derivative,
%  or linear differential operators to be applied before evaluation.
%  Their defaults are 0, meaning that the function itself is evaluated.
%  In the resulting array, the first dimension corresponds to s and
%  the second to t.  

%  last modified 20 July 2006

%  allow different order for 1st three arguments

if nargin >= 3 && isnumeric(bifd) && isa_Lfd(tevalarg)
   temp = bifd;
   bifd = tevalarg;
   tevalarg = sevalarg;
   sevalarg = temp;
end

  if nargin < 5
    tLfd = 0;
  end
  if nargin < 4
    sLfd = 0;
  end

  if ~isa_bifd(bifd)
    error('Argument BIFD is not a bivariate functional data object.');
  end
  ns = length(sevalarg);
  if ~isa_Lfd(sLfd)
    error (['Argument SLfd is neither a functional data object', ...
             ' nor an integer.']);
  end

  if ~isa_Lfd(tLfd)
    error (['Argument TLfd is neither a functional data object', ...
             ' nor an integer.']);
  end

  if strcmp(class(sLfd), 'double')
    if length(sLfd) == 1
      snderiv = round(sLfd);
      if snderiv ~= sLfd
        error('Order of derivative must be an integer');
      end
      if snderiv < 0
        error('Order of derivative must be 0 or positive');
      end
    else
      error('Order of derivative must be a single number');
    end
  else
    sderivcoef  = getcoef(sLfd);
    sderivcoefd = size(sderivcoef);
    snderiv     = sderivcoefd(2);
  end

  if snderiv > 0 && ~strcmp(class(sLfd), 'double')
    derivwtmat = eval(sLfd, sevalarg);
    onerow <- ones(1,snbasis);
    for j = 1:snderiv
      if any(abs(derivwtmat(:,j))) > 1e-7
        sbasismat <- sbasismat +  ...
          (derivwtfdmat(:,j)*onerow) .* ...
          getbasismatrix(sevalarg, sbasisobj, j-1);
      end
    end
  end

  nt = length(tevalarg);

  if strcmp(class(tLfd), 'double')
    if length(tLfd) == 1
      tnderiv = round(tLfd);
      if tnderiv ~= tLfd
        error('Order of derivative must be an integer');
      end
      if tnderiv < 0
        error('Order of derivative must be 0 or positive');
      end
    else
      error('Order of derivative must be a single number');
    end
  else
    tderivcoef  = getcoef(tLfd);
    tderivcoefd = size(tderivcoef);
    tnderiv     = tderivcoefd(2);
  end

  if tnderiv > 0 && ~strcmp(class(tLfd), 'double')
    derivwtmat = eval(tLfd, tevalarg);
    onerow <- ones(1,tnbasis);
    for j = 1:tnderiv
      if any(abs(derivwtmat(:,j))) > 1e-7
        tbasismat <- tbasismat +  ...
          (derivwtfdmat(:,j)*onerow) .* ...
          getbasismatrix(tevalarg, tbasisobj, j-1);
      end
    end
  end

  coef  = bifd.coef;
  coefd = size(coef);
  ndim  = length(coefd);

  switch ndim
  case 2
    evalarray = sbasismat * coef * tbasismat';
  case 3
    ncurves  = coefd(3);
    evalarray = zeros(ns,nt,ncurves);
    for i= 1:ncurves
      evalarray(:,:,i) = sbasismat * coef(:,:,i) * tbasismat';
    end
  case 4
    ncurves  = coefd(3);
    nvar     = coefd(4);
    evalarray = zeros(ns,nt,ncurves,nvar);
    for i = 1:ncurves
      for j = 1:nvar
        evalarray(:,:,i,j) = sbasismat * coef(:,:,i,j) * tbasismat';
      end
    end
  otherwise
    error('coefficient array of improper dimension');
  end
