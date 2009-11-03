function  pcarotstr = varmx_pca(pcastr, nharm, nx)
%  VARMX_PCA  Apply varimax rotation to harmonics in a pca.fd object.
%  Arguments:
%  PCASTR ... Struct object returned by PCA_FD.
%  NHARM  ... Number of harmonics to rotate.
%  NX     ... Number of equally spaced points at which harmonics are
%             evaluated.  Default 501.
%  Return:
%  PCAROTSTR ... PCASTR object containing rotated harmonics.

%  Last modified 30 June 2010 by Jim

  if nargin < 3
    nx = 501;
  end
  
  harmfd    = pcastr.harmfd;
  harmcoef  = getcoef(harmfd);
  coefd     = size(harmcoef);
  ndim      = length(coefd);
  
  varprop0  = pcastr.varprop;
 
  scoresd   = size(pcastr.harmscr);
  
  if nargin < 2
    nharm = scoresd(2);
  end
  
  if nharm > scoresd(2)
    nharm = scoresd(2);
  end

  basisobj = getbasis(harmfd);
  rangex   = getbasisrange(basisobj);
  x        = linspace(rangex(1), rangex(2), nx);
  harmmat  = eval_fd(x, harmfd);
  
  %  If fdmat is a 3-D array, stack into a matrix
  
  if ndim == 3
    harmmatd = size(harmmat);
    harmmat  = permute(harmmat, [1, 3, 2]);
    harmmat  = reshape(harmmat, harmmatd(1) * harmmatd(3), harmmatd(2));
  end
  
  %  compute rotation matrix for varimax rotation of harmmat
  
  rotmat = varmx(harmmat(:,1:nharm));
  
  %  rotate harmonic coefficients 
  
  if ndim == 2
    harmcoef = harmcoef(:,1:nharm);
    harmcoef = harmcoef * rotmat;
  else
    harmcoef = harmcoef(:,1:nharm,:);
    for j = 1:coefd(3)
      harmcoef(:,1:nharm,j) = harmcoef(:,1:nharm,j) * rotmat;
    end
  end
  
  %  modify harmfd object
  
  harmfd = putcoef(harmfd, harmcoef);
  
  %  rotate principal component scores
  
  harmscr = pcastr.harmscr;
  harmscr(:,1:nharm) = harmscr(:,1:nharm) * rotmat;

  %  Compute rotated harmonic variances
  
  propvar = squeeze(sum(harmscr.^2));
  if ndim == 2
      propvar = propvar./sum(propvar);
      propvar = propvar.*sum(varprop0);
  else
      for j=1:coefd(3)
          propvar(:,j) = propvar(:,j)./sum(propvar(:,j));
          propvar(:,j) = propvar(:,j).*sum(varprop0);
      end
  end
  
  %  set up the output pca struct object
  
  pcarotstr         = pcastr;
  pcarotstr.harmfd  = harmfd;
  pcarotstr.harmscr = harmscr;
  pcarotstr.varprop = propvar;
  pcarotstr.values  = pcastr.values;
  pcarotstr.meanfd  = pcastr.meanfd;

