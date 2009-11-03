function [casenames, varnames] = getfdlabels(fdnames) 

%  fdnames must be a cell object of length 3

if iscell(fdnames) && length(fdnames) == 3
    
    %  set up casenames:  fdnames{2} must be a cell object of length 2
    
    if iscell(fdnames{2}) && length(fdnames{2}) == 2
        casenames = fdnames{2};
    else
        %  if not, casenames is empty
        casenames = [];
    end
    
    %  set up varnames:  fdnames{3} must be a cell object of length 2
    
    if iscell(fdnames{3}) && length(fdnames{3}) == 2
        varnames = fdnames{3};
    else
        %  if not, varnames is empty
        varnames = [];
    end
    
else
    
    %  if not, case and variable names are empty
    
    casenames = [];
    varnames  = [];
    
end
    


