function   plot_cca(ccastr, overplt, jcan, flip)
%  PLOT_CCA  Plot a functional canonical correlation analysis object.
%
%  If overplt=1  then each pair of weight functions is plotted in
%     a single plot.  The line types and colours of the
%     'x' and 'y' curves respectively are specified as in plot.str   .
%  If overplt=0  then the weight functions are plotted in separate
%     plots, side by side if a command like par(mfrow=c(2,2)) is
%       used.
%
%  If jcan=0, then all the pairs of variates are plotted.  Otherwise
%     only the variates jcan are plotted (eg if jcan=1, only the leading
%     variate is plotted, if jcan=c(1,3) only the first and third.)
%
%  If flip(j) is 1 then the jth pair of weight functions are multiplied
%     by -1.  If flip is a scalar it is replicated to the necessary length.
%
%  Other arguments are passed to plot_str
%

%  Last modified 3 July 1998

  if nargin < 4
    flip = 0;
  end

  if nargin < 3
    jcan = 0;
  end

  if nargin < 2
    overplt = 0;
  end

  wtfdobj = ccastr.wtfdobj;
  wtcoef = getcoef(wtfdobj);
  if jcan(1) ~= 0
    wtcoef = wtcoef(:,jcan,:);
    wtcoef = permute((-1)^flip * permute(wtcoef, [3, 2, 1]), [3, 2, 1]);
  end
  if overplt
    wtfdobj = putcoef(wtfdobj, wtcoef);
    plot(wtfdobj);
  else
    wtcoef = permute(wtcoef, [1, 3, 2]);
    wtcoefd = size(wtcoef);
    ncan = wtcoefd(2);
    for jj = 1:ncan
      wtfdobjtemp = wtfdobj;
      wtcoefjj = squeeze(wtcoef(:,jj,:,:));
      wtfdobjtemp = putcoef(wtfdobjtemp, wtcoefjj);
      plot(wtfdobjtemp)
    end
  end
