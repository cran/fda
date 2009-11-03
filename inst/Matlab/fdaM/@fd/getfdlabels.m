function [xlabel, ylabel, casenames, varnames] = ...
                  getfdlabels(fdnames, nrep, nvar) 

%  Extract plot labels and, if available, names for each replicate and
%  each variable

%  check fdnames, which must be a list object of length 3

if ~iscell(fdnames)
    error('Argument fdnames is not a list object.');
end
    

if length(fdnames) ~= 3
    error('Argument fdnames is not of length 3.');
end

%  xlabel is fdnames{1} if it has length 1 and is not null
%  otherwise xlabel is names(fdnames)[1]

xlabel = fdnames{1};
%   if length(xlabel) > 1 || isempty(xlabel) xlabel = names(fdnames)[1]
%   if !is.character(xlabel)) xlabel = ''

%  ylabel is fdnames{3} if it has length not equal to nvar and is not null
%  otherwise ylabel is names(fdnames)[3]

ylabel = fdnames{3};
%   if ( (nvar > 1 && length(ylabel) == nvar) || 
%       is.null(ylabel)) ylabel = names(fdnames)[3]
if length(ylabel) > 1
    if ischar(ylabel)
        ylabel = ylabel(1);
    else 
      if iscell(ylabel)
          ylabel = ylabel{1};
      else
          ylabel = '';
      end
    end
end
if ~ischar(ylabel) 
    ylabel = '';
end

%  set up casenames
 
if length(fdnames{2}) == nrep 
    casenames = fdnames{2};
else                             
    casenames = [];
end

%  set up varnames

if length(fdnames{3}) == nvar 
    varnames = fdnames{3};
else                           
    varnames = [];
end


