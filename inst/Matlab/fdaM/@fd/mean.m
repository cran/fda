function meanfd = mean(fdobj)
%  MEAN  Compute mean function for functional observations
%  Argument:
%  FDOBJ ... a functional data object for the observations
%  Return:
%  MEANFD ... a functional data object for the mean

%  last modified 4 March 2009

  if ~isa_fd(fdobj)
    error ('Argument FD is not a functional data object.');
  end

  coef     = getcoef(fdobj);
  coefd    = size(coef);
  ndim     = length(coefd);
  basisobj = getbasis(fdobj);
  if ndim == 2
    coefmean = mean(coef,2);
  else
    nvar = coefd(3);
    coefmean = zeros(coefd(1),1,nvar);
    for j = 1:nvar
      coefmean(:,1,j) = mean(coef(:,:,j),2);
    end
  end

  meanfdnames = getnames(fdobj);
  meanfdnames{2} = 'Mean';
  if iscell(meanfdnames{3})
      meanfdnames{3}{1} = ['Mean ', meanfdnames{3}{1}];
   else
     meanfdnames{3} = ['Mean ', meanfdnames{3}];      
  end

  meanfd.coef     = coefmean;
  meanfd.basisobj = basisobj;
  meanfd.fdnames  = meanfdnames;

  meanfd = class(meanfd, 'fd');

