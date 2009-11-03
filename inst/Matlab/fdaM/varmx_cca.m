function  [ccawtfd1, ccawtfd2, ccavar1, ccavar2, rotmat] = ...
    varmx_cca(ccawtfd1, ccawtfd2, ccavar1, ccavar2, nx)
%  VARMX_CCA  Apply varimax rotation to CCA weight functions
%             and scores
%  Arguments:
%  CCAWTFD1 ... A functional parameter object for the canonical weight
%                functions for the first  set of functions.
%  CCAWTFD2 ... A functional parameter object for the canonical weight
%                functions for the second set of functions.
%  CCAVAR1  ... Canonical variate scores for first  set of functions.
%  CCAVAR2  ... Canonical variate scores for second set of functions.
%  Return:  Arguments after rotation.

%  last modified 31 March 2004

  if nargin < 5
    nx = 50;
  end
  
  wtcoef1 = getcoef(ccawtfd1);
  wtcoef2 = getcoef(ccawtfd2);
  
  basisobj = getbasis(ccawtfd1);
  rangex   = getbasisrange(basisobj);
  x        = linspace(rangex(1), rangex(2), nx);
  ccawtmat1 = eval_fd(x, ccawtfd1);
  ccawtmat2 = eval_fd(x, ccawtfd2);
  ccawtmat  = [ccawtmat1; ccawtmat2];
  %  If fdmat is a 3-D array, stack into a matrix
  %  compute rotation matrix for varimax rotation of harmmat
  rotmat = varmx(ccawtmat);
  %  rotate coefficients and scores
  wtcoef1 = wtcoef1 * rotmat;
  wtcoef2 = wtcoef2 * rotmat;
  %  rotate ccawt objects
  ccawtfd1 = putcoef(ccawtfd1, wtcoef1);
  ccawtfd2 = putcoef(ccawtfd2, wtcoef2);
  %  rotate cca scores
  ccavar1 = ccavar1 * rotmat;
  ccavar2 = ccavar2 * rotmat;


