function ccastr = cca(fdobj, ncan, lambda, Lfd)
%  CCA_FD   Functional canonical correlation analysis with regularization.
%
%  Arguments:
%  FDOBJ    ... Functional data object.  It is assumed that there are
%                two functions for each replication.
%  NCAN     ... Number of pairs of canonical variates to be found. Default 2
%  LAMBDA   ... Smoothing or regularization parameter.
%               If lambda is a 2-vector then the first
%                component will be applied to the 'x' and the second to
%                the 'y' functions.  Default 0.00025
%  LSTR     ... The order of the derivative to be penalized if an integer, or
%                a linear differential operator if a functional data object.
%                Default 2.
%
%  Returns:  A struct object CCASTR with the fields
%  WTFDSTR  ... A functional data object for the canonical
%                        variate weight functions
%  CORR     ... The corresponding set of canonical correlations.
%  VARS     ... An array of the values taken by the canonical variates (ie
%                     the scores on the variates.)  This is a 3-way array
%                     with first dimension corresponding to replicates,
%                     second to the different variates (dimension NCAN)
%                     and third (dimension 2) to the 'x' and 'y' scores.
%

%  Last modified on:  20 July 2006

  if nargin < 4
      Lfd = int2Lfd(2);
  end
  if nargin < 3
      lambda = 0.00025;
  end
  if nargin < 2
      ncan = 2;
  end
  
  %  Check arguments
  
  if ~isa_fd(fdobj)
    error ('First argument is not a functional data object');
  end

  if ~isa_Lfd(Lfd)
    error ('Fourth Argument is not a linear differential operator object.');
  end

  %  Center functions if needed
  
  fdobj = center(fdobj);
  coef  = getcoef(fdobj);
  coefd = size(coef);
  ncoef = length(coefd);
  if ncoef < 3 || coefd(3) == 1
    error('CCA only possible with bivariate functions');
  end
  if coefd(3) > 2
    error('Multiple functions are not permitted.');
  end
  lambda = [lambda, lambda];
  if lambda(1) <= 0 || lambda(2) <= 0
    error('Smoothing parameters must be strictly positive');
  end

  fdbasis   = getbasis(fdobj);
  nbasis    = getnbasis(fdbasis);
  nrep      = coefd(2);
  if nrep < 2
    error('CCA not possible without replications.');
  end

  %   Set up cross product matrices

  Jmat  = eval_penalty(fdbasis, 0);
  Jx    = (Jmat * coef(:,:,1))';
  Jy    = (Jmat * coef(:,:,2))';
  PVxx  = Jx' * Jx./nrep; 
  PVyy  = Jy' * Jy./nrep;
  if any(lambda > 0)
    Kmat  = eval_penalty(fdbasis, Lfd);
    if lambda(1) > 0 
      PVxx  = PVxx + lambda(1) * Kmat;
    end
    if lambda(2) > 0 
      PVyy  = PVyy + lambda(2) * Kmat;
    end
  end
  Vxy   = Jx' * Jy./nrep;
  %  do eigenanalysis
  geigstr = geigen(Vxy, PVxx, PVyy);
  %  set up canonical correlations and coefficients for weight functions
  canwtcoef = zeros([nbasis,ncan,2]);
  canwtcoef(:,:,1) = geigstr.Lmat(:,1:ncan);
  canwtcoef(:,:,2) = geigstr.Mmat(:,1:ncan);
  corrs = diag(geigstr.values);
  corrs = corrs(1:ncan);

  %   Normalize the weight functions

  for j = 1:ncan
    for k = 1:2
      temp = squeeze(canwtcoef(:,j,k));
      temp = temp./sqrt(sum(temp.^2));
      canwtcoef(:,j,k) = temp;
    end
  end

  %  set up final results

  fdnames = getnames(fdobj);
  wtfd    = fd(canwtcoef, fdbasis, fdnames);

  ccavars(:,:,1) = Jx * canwtcoef(:,:,1);
  ccavars(:,:,2) = Jy * canwtcoef(:,:,2);

  ccastr.wtfdobj = wtfd;
  ccastr.corr    = corrs;
  ccastr.vars    = ccavars;
