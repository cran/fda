function varbifd = var(xfd, yfd)
%  VAR  Compute variance and covariance functions .
%  VAR  is normally called with one argument.  In this case,
%     if there are multiple functions, there is a variance-covariate
%     bivariate function computed for each pair of functions.
%     If VAR is called with two arguments, both must be univariate
%     functions, and in this case only the cross-covariance bivariate
%     function is computed.
%  Arguments:
%  XFD ... a functional data object
%  YFD ... an optional second functional data object
%             It is assumed that each object has the same number of
%             replications.
%  Return:
%  VARBIFD  ...  a bivariate functional data object.
%     If VAR is called with only one argument:
%           If the functions are univariate, then there is
%             a single bivariate variance-covariance function.
%           If the functions are multivariate, there is a
%             variance-covariance function for each
%             pair of functions.
%     If VAR is called with two arguments:
%           Both functions must be univariate, and in this case a
%           single cross-covariance function is computed.
%     In either case, the bivariate variance-covariance function(s)
%           have a single replication, so that the third index of
%           the coefficient array is 1.

%  last modified 20 July 2006

  coefx     = getcoef(xfd);
  coefdx    = size(coefx);
  ndimx     = length(coefdx);
  xbasisobj = getbasis(xfd);
  xnbasis   = getnbasis(xbasisobj);
  nrepx     = coefdx(2);
  onen      = ones(nrepx,1);

  fdnames = getnames(xfd);

  if nargin < 2

  %  ----------------------------------------------------------------------
  %  only one argument ... get variance-covariance bivariate function
  %  ----------------------------------------------------------------------

    if ~isa_fd(xfd)
      error ('Argument XFD is not a functional data object.');
    end

    if ndimx == 2
      %  In this case the coefficient array has two dimensions
      temp    = coefx';
      temp    = temp - onen * mean(temp);
      coefvar = (temp' * temp)./nrepx;
    else
     %  In this case the coefficient array has four dimensions,
     %     the third is 1, and the fourth is equal to the number
     %     of pairs of functions.
      nvar  = coefdx(3);
      npair = nvar*(nvar+1)/2;
      coefvar = zeros(xnbasis,xnbasis,1,npair);
      m = 0;
      for i = 1:nvar
        tempi = coefx(:,:,i)';
        tempi = tempi - ones(nrepx,1) * mean(tempi);
        for j = 1:i
          m = m + 1;
          tempj = coefx(:,:,j)';
          tempj = tempj - ones(nrepx,1) * mean(tempj);
          varmat = (tempi' * tempj)./nrepx;
          coefvar(:,:,1,m) = varmat;
        end
      end
    end

    varbifd = bifd(coefvar, xbasisobj, xbasisobj, fdnames);

    return;

  else

  %  ----------------------------------------------------------------------
  %  two arguments ... get covariance bivariate function, if both
  %    FD's are univariate.  Multivariate case not yet set up.
  %  ----------------------------------------------------------------------

    if ~isa_fd(xfd)
      error ('Argument XFD is not a functional data object.');
    end
    if ~isa_fd(yfd)
      error ('Argument YFD is not a functional data object.');
    end

    coefy     = getcoef(yfd);
    coefdy    = size(coefy);
    ndimy     = length(coefdy);
    ybasisobj = getbasis(yfd);

    if coefdx(2) ~= coefdy(2)
      error('Number of replications are not equal.');
    end

    if ndimx == 2 && ndimy == 2

      tempi = coefx';
      tempi = tempi - onen * mean(tempi);
      tempj = coefy';
      tempj = tempj - onen * mean(tempj);
      coefvar = (tempi' * tempj)./nrepx;

      varbifd = bifd(coefvar, xbasisobj, ybasisobj, fdnames);

    else
      error('Both arguments must be univariate.');
    end
  end

