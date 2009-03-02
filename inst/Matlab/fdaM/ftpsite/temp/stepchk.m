function [aval,ind,limwrd] = stepchk(aval, cvec, deltac, limwrd, ind, ...
                                     climit, active, dbgwrd)
%STEPCHK checks the step along a line for producing parameters within the
%  limits specified by BOT and TOP
%  LIMWRD  ...  Logical variable permitting detection that parameter
%               was on the boundary two steps = a row
T   = 1;
F   = 0;
if nargin < 8
  dbgwrd = F;
end
npar = length(deltac);
if nargin < 7
  climit = [-50.*ones(1,npar); 50.*ones(1,npar)];
end
bot    = climit(1,:)';
top    = climit(2,:)';
if aval < 1e-7
  ind = 1;
  return;
end
%  ensure that step does not go beyond lower limit on parameters
stepi   = aval.*deltac;
stepmin = min(stepi);
if any(stepi(active) < bot(active)-cvec(active))
  index   = active(find(stepi(active) < bot(active)-cvec(active) & deltac(active) ~= 0));
  stepnew = min((bot(index)-cvec(index))./deltac(index));
  if dbgwrd
    fprintf('Lower limit reached, old step, new step: %10.4f, %10.4f\n', aval, stepnew);
  end
  aval = stepnew;
  %  check whether lower limit has been reached twice in a row
  if limwrd(1)
    ind = 1;
    return;
  else
    limwrd(1) = T;
  end
else
  limwrd(1) = F;
end
if aval < 1e-7
  ind = 1;
  return;
end
%  ensure that step does not go beyond upper limit on parameters
stepi   = aval*deltac;
stepmax = max(stepi);
if any(stepi(active) > top(active)-cvec(active))
  index   = active(find(stepi(active) > top(active)-cvec(active) & deltac(active) ~= 0));
  stepnew = min((top(index)-cvec(index))./deltac(index));
  if dbgwrd
    fprintf('Upper limit reached, old step, new step: %10.4f, %10.4f\n', aval, stepnew);
  end
  aval = stepnew;
  %  check whether upper limit has been reached twice in a row
  if limwrd(2)
    ind = 1;
    return;
  else
    limwrd(2) = T;
  end
else
  limwrd(2) = F;
end

